## RIMA tutorial:  
Please find our detail tutorial at https://liulab-dfci.github.io/RIMA/.

## Note
if you want to run cibersort in immune infitration module, you should download the source code **CIBERSORT.R** and signature matrix **LM22.txt** from [CIBERSORT website](https://cibersort.stanford.edu).

Make sure you should download **CIBERSORT.R** into RIMA_pipeline/src/immune_infiltration/ folder;\
**LM22.txt** shoube be in RIMA_pipeline/static/cibersort folder.

RIMA has required an overhaul in a lot of the code and requirement of conda environment, we provided the the required dependencies of environment at static/environment/. There could be some bugs due to the specificity of computing clusters used to run this pipeline. We would greatly appreciate if you post any issues to the issues page and will do our best to fix them as quickly as possible. We currently silience the Gene Set Enrichment Analysis (GSEA) function is currently unavaible due to conda environment issue.

## Support 
Create a [GitHub issue](https://github.com/liulab-dfci/RIMA_pipeline/issues).

Or contact main developers:\
Yang Liu: yangliu@ds.dfci.harvard.edu\
Lin Yang: lyang2021nuli@gmail.com 

## Citing RIMA  
If you use RIMA pipeline in your work, please cite [Yang et al.2023](https://www.nature.com/articles/s41596-023-00841-8)
