# Pipeline to get neutral regions to build neutral SFS for vaquita

## 1. GetHiQualCoords_20200602.sh (uses the script obtain_high_qual_coordinates.py)

## 2. exon_distance_calc.sh

## 3. Extract_noCpG_noRepetitive.sh

## 4. get_Coord_file.sh (uses the script obtain_noCpG_noRepetitive_coordinates.py)

## 5. identify_Conserved_Regions.sh

## 6. neutral_Sites2vcf.sh (uses the script obtain_noCpG_noRepetitive_coordinates.py)
