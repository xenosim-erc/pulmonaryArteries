#!/usr/bin/env bash

# Do not change the caller's error-handling behaviour when this file is sourced.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -e
fi

check_openfoam_environment()
{
if ! command -v foamDictionary >/dev/null 2>&1; then
    echo "ERROR: OpenFOAM is not loaded. Load the OpenFOAM v2412 environment and try again."
    return 1
fi

if [[ "${WM_PROJECT_VERSION:-}" != "v2412" ]]; then
    echo "ERROR: This script requires OpenFOAM v2412 (found '${WM_PROJECT_VERSION:-unknown}')."
    return 1
fi
}

require_positive_number()
{
local value="$1"
local description="$2"

if ! awk -v value="$value" 'BEGIN {
    number = "^[+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$"
    exit !(value ~ number && (value + 0) > 0)
}'; then
    echo "ERROR: $description must be a positive number (received '$value')."
    return 1
fi
}

clean_case()
{
read -r -p "Would you like to clean the case? (y/n): " CLEAN_CASE
if [[ "${CLEAN_CASE,,}" == "y" || "${CLEAN_CASE,,}" == "yes" ]]; then
    ./Allclean
fi
}

prepare_case_and_mesh()
{
local template_dir="${MURRAY_TEMPLATE_DIR:-${WM_PROJECT_USER_DIR:-}/pulmScript}"

if [[ ! -d "$template_dir" ]]; then
    echo "ERROR: Cannot find the case template directory: $template_dir"
    echo "Set MURRAY_TEMPLATE_DIR to the directory containing system, constant, 0, run.slurm and Allrun."
    return 1
fi

read -p "Enter geometry file name:" GEOM
echo $GEOM

cp -r \
    "$template_dir/system" \
    "$template_dir/constant" \
    "$template_dir/0" \
    "$template_dir/run.slurm" \
    "$template_dir/Allrun" \
    .

surfaceFeatureEdges $GEOM.stl $GEOM.fms

##Creating Mesh
read -p "Enter largest cell size:" cMax
read -p "Enter smallest cell size:" cMin
require_positive_number "$cMax" "Largest cell size" || return 1
require_positive_number "$cMin" "Smallest cell size" || return 1
#read -p "Enter number of boundary layers:" nBound
#read -p "Enter boundary cell size:" cBound

sed -i "s/maxCellSize.*/maxCellSize $cMax;/" system/meshDict
sed -i "s/minCellSize.*/minCellSize $cMin;/" system/meshDict

sed -i "s/surfaceFile.*/surfaceFile $GEOM.fms;/" system/meshDict

echo meshing...
cartesianMesh >>log.cartesianMesh
echo checking mesh...
checkMesh >>log.checkMesh

if grep -q "^Mesh OK\." log.checkMesh; then
    echo Mesh is a work of art - Behold ye mighty and despair
    echo Number of cells is:
    grep "cells:" log.checkMesh 
else
    echo Mesh is a disgusting failure
    grep "^Failed" log.checkMesh
    echo Try other cell sizes
fi

read -p "are you happy with the mesh? (y/n)" mOK
if [[ $mOK == "n" ]]; then
    return 1
fi

    
echo creating patches...
autoPatch 50 -overwrite >>log.autoPatch

auto_patch_count=$(awk '
    /^[[:space:]]*auto[0-9]+[[:space:]]*$/ { count++ }
    END { print count + 0 }
' constant/polyMesh/boundary)
outlet_count=$((auto_patch_count - 2))

echo "Number of outlets found: $outlet_count"

read -r -p "Is this number of outlets correct? (y/n): " outlets_correct
if [[ "${outlets_correct,,}" != "y" && "${outlets_correct,,}" != "yes" ]]; then
    echo "The outlet count is incorrect. Please clean up the geometry and try again."
    return 1
fi



##Write patches to file in size order
patches=($(awk '
/auto[0-9]+/ {name=$1}
/nFaces/ {gsub(";", "", $2); print name, $2}
' constant/polyMesh/boundary |
sort -k2 -nr |
awk '{print $1}'))

echo "${patches[0]}"   # largest
echo "${patches[1]}"   # second largest

largest="${patches[0]}"   # largest
inlet="${patches[1]}"   # second largest

echo $largest

#largest patch must be wall
#second largest must be inlet
#safe assumption for pulmonary geometry

sed -i "s/\<$largest\>/wall/" constant/polyMesh/boundary
createPatch -overwrite

sed -i "s/\<$inlet\>/inlet/" constant/polyMesh/boundary

read -p "Enter scale factor to put mesh in metres:" scale 
require_positive_number "$scale" "Scale factor" || return 1
transformPoints -scale "$scale"
}


##Code reads in all of the patches and then writes a function object to output all areas to postProcessing. These are then read back in

apply_murray_rule()
(
CASE_DIR="${1:-.}"
FIELD="${2:-U}"

cd "$CASE_DIR" || return 1

BOUNDARY_FILE="constant/polyMesh/boundary"
if [ ! -f "$BOUNDARY_FILE" ]; then
    echo "ERROR: Cannot find $BOUNDARY_FILE"
    echo "Run this script from (or point it at) the top level of an OpenFOAM case."
    return 1
fi

echo "Reading patch names from $BOUNDARY_FILE ..."

# Try foamDictionary first (most reliable), fall back to an awk parse
PATCHES=$(foamDictionary -entry boundary -keywords "$BOUNDARY_FILE" 2>/dev/null || true)

if [ -z "$PATCHES" ]; then
    PATCHES=$(awk '
        /^[[:space:]]*[A-Za-z0-9_.\-]+[[:space:]]*$/ {
            gsub(/^[ \t]+|[ \t]+$/, "", $0)
            candidate=$0
        }
        /nFaces/ { print candidate }
    ' "$BOUNDARY_FILE")
fi

if [ -z "$PATCHES" ]; then
    echo "ERROR: Could not parse any patch names from $BOUNDARY_FILE"
    return 1
fi

echo "Found patches:"
echo "$PATCHES" | sed 's/^/  - /'
echo

# Build a function object dictionary, one surfaceFieldValue per patch
mkdir -p system
FUNC_DICT="system/patchAreasDict"

{
cat <<'HEADER'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      patchAreasDict;
}
HEADER
echo "functions"
echo "{"
for p in $PATCHES; do
cat <<EOF
    area_${p}
    {
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        enabled         true;
        writeControl    writeTime;
        writeFields     false;
        surfaceFormat   none;
        regionType      patch;
        name            ${p};
        operation       areaIntegrate;
        writeArea       true;
        fields          (U);
    }

EOF
done
echo "}"
} > "$FUNC_DICT"

echo "Running postProcess -dict $FUNC_DICT -latestTime ..."
postProcess -dict "$FUNC_DICT" -latestTime > patchAreas.log 2>&1

echo
echo "==================== Patch areas ===================="
printf "%-25s %s\n" "Patch" "Area [m^2]"
printf "%-25s %s\n" "-----" "----------"

for p in $PATCHES; do
    f=$(find "postProcessing/area_${p}" -name "surfaceFieldValue.dat" 2>/dev/null | sort | tail -1)
    if [ -f "$f" ]; then
        # Last data line, second column = writeArea output (Area)
        area=$(tail -n 1 "$f" | awk '{print $2}')
        printf "%-25s %s\n" "$p" "$area"
    else
        printf "%-25s %s\n" "$p" "N/A (check patchAreas.log)"
    fi
done
echo "======================================================="
echo
echo "Full solver output logged to: patchAreas.log"
echo "Raw per-patch data in:        postProcessing/area_<patchName>/<time>/surfaceFieldValue.dat"


###Code from here reads in the outlet patch areas -> converts to decimal -> puts area to 3/2 -> sums these exponentiated area -> additionally writes each patches exponentiated area to a variable 

sumS=0
outlet_area_count=0

for p in $PATCHES; do
    if [[ "$p" == auto* ]];then
	echo "Processing $p" 
    f=$(find "postProcessing/area_${p}" -name "surfaceFieldValue.dat" 2>/dev/null | sort | tail -1)
    if [ ! -f "$f" ]; then
        echo "ERROR: No area data was found for outlet patch '$p'. Check patchAreas.log."
        return 1
    fi

    # Last data line, second column = writeArea output (Area)
    area=$(tail -n 1 "$f" | awk '{print $2}')
	require_positive_number "$area" "Area for outlet patch '$p'" || return 1
	
    areaDec=$(echo "$area" | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}')
	require_positive_number "$areaDec" "Converted area for outlet patch '$p'" || return 1
    declare -g "area_$p=$areaDec"
    echo The area of this patch is: $areaDec
    
    aExp=$(echo "$areaDec" | bc -l)
    declare -g "areaS_$p=$aExp"
    
    echo "The area raised to the power 3/2 is: $aExp"
    sumS=$(echo "$sumS + $aExp" | bc -l)
    outlet_area_count=$((outlet_area_count + 1))
    
    fi
done

if (( outlet_area_count == 0 )); then
    echo "ERROR: No outlet patches were found for the Murray-rule calculation."
    return 1
fi

require_positive_number "$sumS" "Sum of outlet areas raised to the power 3/2" || return 1

echo $sumS

read -p "Enter total distal resistance (typically 20000000 for human):" Rd
read -p "Enter total proximal resistance (typically 3300000 for human):" Rp
read -p "Enter total downstream capacitance (typically 0.000000035 for human):" C
require_positive_number "$Rd" "Total distal resistance" || return 1
require_positive_number "$Rp" "Total proximal resistance" || return 1
require_positive_number "$C" "Total downstream capacitance" || return 1



##Writing p file

cat > 0/p <<HEADER
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
HEADER

##Applying Area Rule
for p in $PATCHES; do
    if [[ "$p" == auto* ]];then
	echo "Processing $p"
	var="areaS_$p"
	AsRatio=$(echo "${!var} / $sumS" | bc -l)
    
	#Capacitances
	C1=$(echo "$AsRatio * $C" | bc -l)
	
	#Proximal Resistances
	RpInv=$(echo "1 / $Rp" | bc -l)
	Rp1Inv=$(echo "$AsRatio * $RpInv" | bc -l)
	Rp1=$(echo "1 / $Rp1Inv" | bc -l)
	
	#Distal Resistances
	RdInv=$(echo "1 / $Rd" | bc -l)
	Rd1Inv=$(echo "$AsRatio * $RdInv" | bc -l)
	Rd1=$(echo "1 / $Rd1Inv" | bc -l)
	
	cat >> 0/p <<EOF
    $p
    {
       type        windkesselPressure;
       R           $Rd1;
       Rch         $Rp1;
       C           $C1;
       rho         1060;
       value       uniform 0;
    }


EOF
    fi
    
    if [[ "$p" == inlet* ]];then
	cat >>0/p <<EOF
       inlet
    {
        type            zeroGradient;
    }


EOF
fi
if [[ "$p" == wall* ]];then
    cat >>0/p <<EOF
      "wall.*"
    {
        type            zeroGradient;
    }

EOF
fi
done
cat >>0/p <<EOF
      processor
    {
      type processor;
    }
}

// ************************************************************************* //
EOF
)

write_velocity_field()
{
read -p "Enter h for human or p for pig:" ANIMAL
read -p "How many heartbeats would you like to run?" HB
require_positive_number "$HB" "Number of heartbeats" || return 1

if [[ "$ANIMAL" == p ]];then
    RT=$(echo "$HB * 0.66" | bc -l)
sed -i "s/^endTime.*/endTime $RT;/" system/controlDict

cat > 0/U <<EOF
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
  // #includeEtc "caseDicts/setConstraintTypes"

  processor
    {
      type processor;
    }
  
    inlet
    {
        type     flowRateInletVelocity;
	volumetricFlowRate  table
(			
(	0	0	)
(	0.01	0.0000031	)
(	0.02	0.000011	)
(	0.03	0.0000262	)
(	0.04	0.0000489	)
(	0.05	0.0000785	)
(	0.06	0.0001114	)
(	0.07	0.0001431	)
(	0.08	0.0001715	)
(	0.09	0.0001944	)
(	0.1	0.0002109	)
(	0.11	0.0002218	)
(	0.12	0.0002283	)
(	0.13	0.0002312	)
(	0.14	0.0002309	)
(	0.15	0.0002288	)
(	0.16	0.0002253	)
(	0.17	0.0002205	)
(	0.18	0.0002143	)
(	0.19	0.0002074	)
(	0.2	0.0001998	)
(	0.21	0.0001916	)
(	0.22	0.0001829	)
(	0.23	0.0001738	)
(	0.24	0.0001643	)
(	0.25	0.0001544	)
(	0.26	0.0001442	)
(	0.27	0.0001335	)
(	0.28	0.0001219	)
(	0.29	0.0001099	)
(	0.3	0.0000975	)
(	0.31	0.0000851	)
(	0.32	0.0000728	)
(	0.33	0.0000615	)
(	0.34	0.0000516	)
(	0.35	0.0000439	)
(	0.36	0.0000381	)
(	0.37	0.0000343	)
(	0.38	0.0000321	)
(	0.39	0.0000313	)
(	0.4	0.0000318	)
(	0.41	0.0000327	)
(	0.42	0.0000335	)
(	0.43	0.0000337	)
(	0.44	0.0000332	)
(	0.45	0.0000323	)
(	0.46	0.000031	)
(	0.47	0.0000297	)
(	0.48	0.0000285	)
(	0.49	0.0000273	)
(	0.5	0.0000263	)
(	0.51	0.0000252	)
(	0.52	0.000024	)
(	0.53	0.0000231	)
(	0.54	0.0000226	)
(	0.55	0.0000226	)
(	0.56	0.000023	)
(	0.57	0.0000241	)
(	0.58	0.0000253	)
(	0.59	0.0000255	)
(	0.6	0.0000246	)
(	0.61	0.0000224	)
(	0.62	0.0000187	)
(	0.63	0.0000143	)
(	0.64	0.0000094	)
(	0.65	0.0000051	)
(	0.66	0	)
			);
	outOfBounds repeat;

	  value uniform (0 0 0);
	
    }

    "wall.*"
      {
	type    noSlip;
      }
    
    "auto.*"
    {
        type            inletOutlet;
	value           uniform (0 0 0 );
	inletValue      uniform (0 0 0);
    }
    

}


// ************************************************************************* //

EOF

    

else
       RT=$(echo "$HB * 0.86" | bc -l)
       sed -i "s/^endTime.*/endTime $RT;/" system/controlDict
    cat > 0/U <<EOF
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
  // #includeEtc "caseDicts/setConstraintTypes"

  processor
    {
      type processor;
    }
  
    inlet
    {
        type     flowRateInletVelocity;
	volumetricFlowRate  table
	  (
(	0	0	        )
(	0.01	0.00001105	)
(	0.02	0.0000271	)
(	0.03	0.00004877	)
(	0.04	0.00007565	)
(	0.05	0.00010639	)
(	0.06	0.00013965	)
(	0.07	0.00017165	)
(	0.08	0.00020055	)
(	0.09	0.00022555	)
(	0.1	0.00024572	)
(	0.11	0.00025968	)
(	0.12	0.00026841	)
(	0.13	0.00027223	)
(	0.14	0.00027223	)
(	0.15	0.00026833	)
(	0.16	0.00026331	)
(	0.17	0.00025648	)
(	0.18	0.00024827	)
(	0.19	0.00023876	)
(	0.2	0.00022835	)
(	0.21	0.00021752	)
(	0.22	0.00020617	)
(	0.23	0.0001955	)
(	0.24	0.00018423	)
(	0.25	0.00017217	)
(	0.26	0.00016029	)
(	0.27	0.00014721	)
(	0.28	0.00013298	)
(	0.29	0.00011854	)
(	0.3	0.00010236	)
(	0.31	0.00008436	)
(	0.32	0.00006639	)
(	0.33	0.00004944	)
(	0.34	0.00003517	)
(	0.35	0.00002115	)
(	0.36	0.00001059	)
(	0.37	0.0000031	)
(	0.38	-0.00000315	)
(	0.39	-0.00000894	)
(	0.4	-0.00001249	)
(	0.41	-0.00001523	)
(	0.42	-0.00001714	)
(	0.43	-0.00001774	)
(	0.44	-0.00001674	)
(	0.45	-0.0000154	)
(	0.46	-0.00001286	)
(	0.47	-0.00001032	)
(	0.48	-0.00000723	)
(	0.49	-0.0000049	)
(	0.5	-0.00000203	)
(	0.51	-0.00000023	)
(	0.52	0.00000115	)
(	0.53	0.00000323	)
(	0.54	0.00000327	)
(	0.55	0.00000327	)
(	0.56	0.00000327	)
(	0.57	0.00000327	)
(	0.58	0.00000327	)
(	0.59	0.00000327	)
(	0.6	0.00000327	)
(	0.61	0.00000327	)
(	0.62	0.00000327	)
(	0.63	0.00000327	)
(	0.64	0.00000327	)
(	0.65	0.00000272	)
(	0.66	0.00000327	)
(	0.67	0.00000268	)
(	0.68	0.00000268	)
(	0.69	0.00000268	)
(	0.7	0.00000268	)
(	0.71	0.00000268	)
(	0.72	0.00000268	)
(	0.73	0.0000021	)
(	0.74	0.0000021	)
(	0.75	0.0000021	)
(	0.76	0.0000021	)
(	0.77	0.0000021	)
(	0.78	0.0000021	)
(	0.79	0.0000021	)
(	0.8	0.0000021	)
(	0.81	0.00000152	)
(	0.82	0.0000002	)
(	0.83	-0.00000117	)
(	0.84	-0.00000314	)
(	0.85	-0.0000049	)
(	0.86	-0.0000049	)
	  
	  );
	outOfBounds repeat;
	  value uniform (0 0 0);
	
    }

    "wall.*"
      {
	type    noSlip;
      }
    
    "auto.*"
    {
        type            inletOutlet;
	value           uniform (0 0 0 );
	inletValue      uniform (0 0 0);
    }
    

}


// ************************************************************************* //
EOF
fi
}

cleanup_postprocessing()
{
rm -rf -- \
    postProcessing \
    log.cartesianMesh \
    log.checkMesh \
    log.autoPatch \
    patchAreas.log
}

main()
(
    local case_dir="${1:-$PWD}"

    check_openfoam_environment

    if [[ ! -d "$case_dir" ]]; then
        echo "ERROR: Case directory does not exist: $case_dir"
        return 1
    fi

    cd "$case_dir" || return 1

    clean_case
    prepare_case_and_mesh
    apply_murray_rule .
    write_velocity_field
    cleanup_postprocessing
)

# Run the complete workflow only when this file is executed directly.
# When sourced, the functions above remain available in the current shell.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
