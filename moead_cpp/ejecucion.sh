#!/bin/bash

#cambiar datos prueba.

cd ficheros
 for j in $(ls);
 do
  echo $j
   ./moea $j 90 0 50 2 1 3000 200
 done


