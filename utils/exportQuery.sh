#!/bin/bash
set -ex

Rscript exportPartitionedTables.R Association/data.csv Study,Dataset,Cancer_Type,Sample Association
Rscript exportPartitionedTables.R Exposure/data.csv Study,Dataset,Cancer_Type,Signature_set_name Exposure
Rscript exportPartitionedTables.R Seqmatrix/data.csv Study,Cancer_Type,Sample,Dataset,Profile Seqmatrix
Rscript exportPartitionedTables.R Signature/data.csv Profile Signature
Rscript exportPartitionedTables.R Signature/data.csv Profile,Signature_set_name Signature

