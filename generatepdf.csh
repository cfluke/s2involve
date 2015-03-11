#!/bin/csh -f 


# s2involve: S2PLOT INteractive VOLumetric Visualisation Environment    
# *************************************************************************
#   Copyright (C) 2015  Christopher Fluke (cfluke@swin.edu.au)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ***************************************************************************/


set template=./template
set publish=./publish

set oname = `echo $1 | tr '_' '-'`
echo Creating 3D-PDF for object: $oname


chmod 755 packagehead.csh
foreach line ("`cat packagehead.csh`")
  set argv = ($line)
  set fname=$1
  set histo=$2
end

set res=`file $fname.tga`
set split = ($res:as/ / /)
set len = $#split

@ hidx = $len
@ widx = $len - 2
@ height = $split[$hidx]

set width = `echo "$split[$widx]/2" | bc`

set w1 = `echo "$width-10" | bc`
set o1 = `echo "$width+5" | bc`

if ($1 == '') then
   echo 'Usage: makemypdf.csh <object name>'

   echo ""

   echo S2PATH: $S2PATH
   echo S2ARCH: $S2ARCH
   echo S2PLOT_PRCDRIVER: $S2PLOT_PRCDRIVER
   echo S2PLOT_IMPATH: $S2PLOT_IMPATH

   echo ""
   echo pdflatex: `which pdflatex`
   echo convert: `which convert`
   echo read: `which read`
   echo cut: `which cut`
   echo awk: `which awk`
   echo sed: `which sed`

   echo ""
   echo template directory: $template
   ls $template

   exit
endif


set dname = DIR-$oname
set pname = $oname\-poster
set cname = $oname\-cbar

rm -irf $dname
mkdir $dname

cp ${template}/template.tex $dname/$oname.tex
cp ${template}/publish.tex $dname/publish-$oname.tex
cp ${template}/s2views.txt  $dname
cp ${template}/s2plot-prc.js $dname
cp s2plotprc.pdf $dname
cp s2direct.prc $dname
cp s2direct.map $dname
cp ${fname} $dname/$oname-meta.txt

set height1=`echo "(${histo}*${height})/1" |bc`
set height2=`echo "${height}-${height1}" |bc`
set cst1 = "${w1}x${height2}+${o1}+0"
set cst2 = "${w1}x${height1}+${o1}+${height2}"
convert -crop ${cst1} ${fname}.tga $dname/$pname.png
convert -crop ${cst2} ${fname}.tga $dname/$cname.png

cd $dname


# Deal with the poster image
set pwidth=`file $pname.png | cut -f2 -d, | awk '{print $1}'`
set pheight=`file $pname.png | cut -f2 -d, | awk '{print $3}'`


if ($pheight > $pwidth) then
   set string="scale=2;$pwidth/$pheight"
   set pratio=`echo $string | bc`
   sed -i "" 's/OBJNAME/'$oname'/' $oname.tex
   sed -i "" 's/HEIGHT/1.0/' $oname.tex
   sed -i "" 's/WIDTH/'$pratio'/' $oname.tex
   sed -i "" 's/OBJNAME/'$oname'/' publish-$oname.tex
   sed -i "" 's/HEIGHT/1.0/' publish-$oname.tex
   sed -i "" 's/WIDTH/'$pratio'/' publish-$oname.tex
else
   set string="scale=2;$pheight/$pwidth"
   set pratio=`echo $string | bc`
   sed -i "" 's/OBJNAME/'$oname'/' $oname.tex
   sed -i "" 's/HEIGHT/'$pratio'/' $oname.tex
   sed -i "" 's/WIDTH/1.0/' $oname.tex
   sed -i "" 's/OBJNAME/'$oname'/' publish-$oname.tex
   sed -i "" 's/HEIGHT/'$pratio'/' publish-$oname.tex
   sed -i "" 's/WIDTH/1.0/' publish-$oname.tex
endif

# Deal with the colour bar/histogram image
set cwidth=`file $cname.png | cut -f2 -d, | awk '{print $1}'`
set cheight=`file $cname.png | cut -f2 -d, | awk '{print $3}'`

if ($cheight > $cwidth) then
   set string="scale=2;$cwidth/$cheight*$pratio"
   set ratio=`echo "$string" | bc`
   sed -i "" 's/CHGHT/'$pratio'/' $oname.tex
   sed -i "" 's/CWDTH/'$ratio'/' $oname.tex
   sed -i "" 's/CHGHT/'$pratio'/' publish-$oname.tex
   sed -i "" 's/CWDTH/'$ratio'/' publish-$oname.tex
else
   set string="scale=2;$cheight/$cwidth*$pratio"
   set ratio=`echo "$string" | bc`
   sed -i "" 's/CHGHT/'$ratio'/' $oname.tex
   sed -i "" 's/CWDTH/'1.0'/' $oname.tex
   sed -i "" 's/CHGHT/'$ratio'/' publish-$oname.tex
   sed -i "" 's/CWDTH/'1.0'/' publish-$oname.tex
endif



pdflatex publish-$oname.tex
pdflatex publish-$oname.tex

open -a "adobe reader" publish-$oname.pdf
open -a "adobe reader" $oname.pdf


