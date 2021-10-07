#!/bin/bash

if [ "$#" = 2 ]
then
    old_name="$1"
    new_name="$2"
    if [ ! -d $old_name ]
    then
	echo "no tool called "$old_name" exists."
    elif [ -d $new_name ]
    then
	echo "there is already a tool called "$new_name"."
    else
	git mv  $old_name $new_name
	git mv $new_name/$old_name.cpp $new_name/$new_name.cpp  
	git mv $new_name/$old_name.h $new_name/$new_name.h
	sed -i "s/"$old_name"/"$new_name"/g" $new_name/*
	sed -i "s/"$old_name"/"$new_name"/g" Unity.cpp
	sed -i "s/"$old_name"/"$new_name"/g" Factory/Factory.cpp
    fi
else
    echo "usage is ./renameTool <old_name> <new_name>"
fi
