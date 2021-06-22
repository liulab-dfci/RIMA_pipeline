#!/bin/bash/



for i in "static/environment"/*.yml
do
  echo "Creating all environments for the pipeline"  
  mamba  env create -f $i
done    
    
    
	 
	 
	 
	 
	    
