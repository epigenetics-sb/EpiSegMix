#!/usr/bin/env bash
while getopts "n:o:i:" OPTION
do
        case $OPTION in
                n)
                        name="$OPTARG" ;;
                o)
                        output_dir=$OPTARG ;;

                i)
                        plot_statistics=$OPTARG ;;
  esac
done

output_file=${output_dir}/${name}.html
touch $output_file
cat << EOF > ${output_file}
${name}
<!DOCTYPE html>
<html>
<head>
<style>
  img {
    max-width: 70%;
    height: auto;
    display: block;
    margin: 0 auto;
  }

  .text-above-figure {
    text-align: center;
    font-weight: bold;
  }
</style>
</head>
<body>
<h1><center>EpiSegMix report</center></h1><br>
<center>
<div class="text-above-figure"><h2>Segmentation</h2><br><a>State colors<br><img src="${name}-stateColors.png"></a></div>
<br>
<div class="text-above-figure"><h2>Average methylation</h2><br><img src="${name}-averageMethylation.png"></div>
<br>
<div class="text-above-figure"><h2>Characteristics</h2><br>
<table>
<tr>
  <td valign="top"><a>Transition matrix<br><img src="${name}-transitionMatrix.png"></a></td>
  <td valign="top"><a>Average length<br><img src="${name}-stateLength.png"></a></td>
  <td valign="top"><a>Coverage<br><img src="${name}-stateMembership.png"></a></td>
</tr>
</table>
</div>
</center>
</body>
</html>
EOF