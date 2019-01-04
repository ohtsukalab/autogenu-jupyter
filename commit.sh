#!/bin/bash
  REMOVE=($(git status | grep deleted | awk '{ print $2}'))
    for i in "${REMOVE[@]}"
    do
      :
      git rm $i
    done

  git add --all .
  git commit -m "$2"
git push "$1" ros-kinetic-devel