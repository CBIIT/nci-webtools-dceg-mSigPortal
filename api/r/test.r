#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
fn <- paste(args[1],"()")

hello <- function(text=args[2]){

    return(text)
}

add <- function(num1=args[2],num2=args[3]){
    return(as.numeric(num1)+as.numeric(num2))
}

bye <- function(){

    return(TRUE)
}

eval(parse(text=fn))
