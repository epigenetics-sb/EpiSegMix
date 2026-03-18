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

if [ "${plot_statistics}" = "False" ] || [ "${plot_statistics}" = "false" ]; then

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
<div class="text-above-figure"><h2>Segmentation</h2><br><a>State colors<br><img src="${name}-state-colors.png"></a></div>
<br>
<div class="text-above-figure"><h2>Normalized counts</h2><br><img src="${name}-normEmission.png"></div>
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

else 

output_file=$(realpath -s ${output_dir}/${name}.html)
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
<div class="text-above-figure"><h2>Segmentation</h2><br><a>State colors<br><img src="${name}-state-colors.png"></a></div>
<br>
<div class="text-above-figure"><h2>Normalized counts</h2><br><img src="${name}-normEmission.png"></div>
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
<br>
<div class="text-above-figure"><h2>Emission distribution</h2><br><img src="${name}-stateDistribution.png"></div>
<br>
<div class="text-above-figure"><h2>Input characteristics</h2><br>
<table>
<tr>
  <td valign="top">Correlation<br><img src="${name}-correlation.png"></td>
  <td valign="top">Distributions<br><img src="${name}-histogram.png"></td>
</tr>
</table>
</div>
</center>
</body>
</html>
EOF

fi 