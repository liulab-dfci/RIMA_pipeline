#!/bin/bash

str2="AWS"
str1="GCP"
while getopts ":p:" opt; do
  case $opt in
    p)
      if [ "$OPTARG" == "$str2" ]; then
          echo "Installing the RIMA environment on: $OPTARG"
	  
	  #install the mamba as conda package builder
	  conda install -c conda-forge mamba
      
	  for i in "static/environment/AWS"/*.yml
          do
            echo "Creating all environments for the pipeline"
            mamba  env create -f $i
          done

	  conda activate RIMA
          sudo apt install libfontconfig1 libxrender1

      elif [ "$OPTARG" == "$str1" ]; then
	  echo "Installing the RIMA environment on: $OPTARG"
	  
	  #install the mamba as conda package builder
          conda install -c conda-forge mamba

          for i in "static/environment/GCP"/*.yml
          do
            echo "Creating all environments for the pipeline"
            mamba  env create -f $i
          done

      else
          echo "$OPTARG is not a support argument"
	  exit 1
      fi
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done
