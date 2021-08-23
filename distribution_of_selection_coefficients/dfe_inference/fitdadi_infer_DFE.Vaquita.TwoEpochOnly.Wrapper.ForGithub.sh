######### DFE wrapper:



#vaquitaDFE inference
conda activate DFEInference # python 2.7 with numpy and scipy 

todaysdate=`date +%Y%m%d`

gitdir=/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/VaquitaDFE/
scriptdir=$gitdir/dfe_inference

resultsdir=/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/

synSFS=$resultsdir/easySFS/synonymous/projection-20200928-hetFilter-0.75/dadi/GoC-24.sfs # syn dir
nonsynSFS=$resultsdir/easySFS/missense/projection-20200928-hetFilter-0.75/dadi/GoC-24.sfs # missense dir

Lcds=30944160 # total CDS sites (not "exon" sites bc those include UTRs) passing filters


mutationRates='2.2e-09,5.8e-09,1.08e-08' # final set of mutation rates for Jacqueline

ns_s_ratios=2.31
exonicMutationRateScalingFactorexonicMutationRateScalingFactor=1 # setting to 1 for these expts 


###########  two epoch ########################
script=fitdadi_infer_DFE.Vaquita.TwoEpoch.GammaOnly.ForGithub.py # this is 2-epoch and gamma dist specific! will redo demographic inference but only for 2epoch model

for NS_S_ScalingFactor in $ns_s_ratios
do
echo -e "beginning with NS:S scaling factor: $NS_S_ScalingFactor"
outprefix=$resultsdir/dfe_inference/two_epoch_model/inference_${todaysdate}_TwoEpochModel_MuGrid_nsScaleFac_${NS_S_ScalingFactor}_FINAL_FOR_MANUSCRIPT # changed so you now input list of mutation rates since dont have to redo demog based on mu rate
 
python $scriptdir/$script \
--syn_input_sfs $synSFS \
--Lcds $Lcds \
--nonsyn_input_sfs $nonsynSFS \
--outprefix $outprefix \
--mutationRateList $mutationRates \
--exonicMutationRateScalingFactorexonicMutationRateScalingFactor $exonicMutationRateScalingFactorexonicMutationRateScalingFactor \
--NS_S_ScalingFactor $NS_S_ScalingFactor 2> ${outprefix}.stdout.log

done


conda deactivate