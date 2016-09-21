#!/bin/bash
#Program:
#    to compare the effect of different a1 and a2
#History:
#2016/9/18 Tingke
#Version: V0.0
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:~/bin
export PATH
g++ *.cpp
for (( i=0;i !=6;++i)) do
	for (( j=0;j !=6;++j)) do
		./a.out $i $j;
		matlab -nodisplay < despersion_relation.m;
		a="ej_"$i"_"$j"_mesh.eps"
		mv ej_mesh.eps $a;
		b="ej_"$i"_"$j"_contour.eps"
		mv ej_contour.eps $b;
		c="ek_"$i"_"$j"_mesh.eps"
		mv ek_mesh.eps $c;
		d="ek_"$i"_"$j"_contour.eps"
		mv ek_contour.eps $d;
	done
done
mv *.eps ./compare_of_a1_a2;
echo -e "The comparation is over"
