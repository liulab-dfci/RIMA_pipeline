#!/bin/bash/



for i in "static/environments"/*.yml
do
  echo "Creating all environments for the pipeline"  
  mamba  env create -f $i
done    
    
    
	 
	 
	 
	 
	    
