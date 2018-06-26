
#This is a short manual how to run ADVOCATE

### ANALYSIS
library(ADVOCATE)

#Load files
infiles = c("ADVOCATE_INPUT_DATA/trainDataADVOCATE.rda",
            "ADVOCATE_INPUT_DATA/input_cumc_bulk_expmat_vst_z.rds")
outfile_prefix =  "Res_ADV_cumcBulk"

bulkexp = readRDS(infiles[2])

## Fractions -------------------------------
## Prediction  of fractions (2 compartments) 
print('frac2')
res.prop = predict_bulk(infiles[1], bulkexp, epsilon = 0.001 )
print(res.prop)

## Prediction of fractions (3 compartments) ####
print('frac3')
res.prop3 = predict_bulk_3comp(infiles[1], bulkexp, epsilon = 0.001 )

### Virtual expressiong ---------------------

## Predict virtual gene expression (2 compartments)
print('ve2')
load(infiles[1])
lcmexp = expmat; rm(expmat)
res.vexp = calCellTypeExpression(lcmexp, deg, fc, pval,sampleInfo, bulkexp, res.prop, method = 'lcm')

## Predict virtual gene expression (3 compartments)
print('ve3')
load(infiles[1])
lcmexp = expmat; rm(expmat)
res.vexp3 = calCellTypeExpression_3comp(lcmexp, deg, fc, pval,sampleInfo,bulkexp, res.prop3, method = 'lcm')

## Output ----------------------------------

wtable = function(...)  write.table(..., sep = "\t", quote = F, col.names = T, row.names =F)

wtable(cbind(rownames(res.prop), res.prop[,1:2]), 
       file = paste0("ADVOCATE_OUTPUT_DATA/",outfile_prefix,"_2c_frac.tsv"))

wtable(cbind(rownames(res.prop3), res.prop3[,1:3]), 
       file = paste0("ADVOCATE_OUTPUT_DATA/",outfile_prefix,"_3c_frac.tsv"))

wtable(cbind(rownames(res.vexp), res.vexp), 
       file = paste0("ADVOCATE_OUTPUT_DATA/",outfile_prefix,"_2c_ve.tsv"))

wtable(cbind(rownames(res.vexp3), res.vexp3), 
       file = paste0("ADVOCATE_OUTPUT_DATA/",outfile_prefix,"_3c_ve.tsv"))

print('END')
