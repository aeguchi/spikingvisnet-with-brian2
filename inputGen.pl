#!/usr/bin/perl
#	shapeGen.pl
#
#	Created by Akihiro Eguchi on 10/30/12.
#
#	Using a PerlMagick library, this 

use Image::Magick;
use Math::Trig;

# setting

$cs = 30;	#canvas size

$img = Image::Magick->new(size=>($cs."x".$cs));
$img->ReadImage('canvas:rgb(255, 255, 255)');
$imgTmp1 = $img->clone();
$size = 8;

system ("mkdir images");
system ("mkdir images/training_gen/");

for($y=$size/2; $y<($cs-($size*1.5)); $y+=1){
	for($x=$size/2; $x<($cs-$size*1.5); $x+=1){
		$img = $imgTmp1->clone();
		#horizontal for L
		$img->Draw(primitive=>'line',stroke=>'black', strokewidth=>'1', fill=>'none', points=>($x).",".($y+$size). " ".($x+$size).",".($y+$size));
		#vertical for L	
		$img->Draw(primitive=>'line',stroke=>'black', strokewidth=>'1', fill=>'none', points=>($x).",".($y). " ".($x).",".($y+$size));

		$img->Write("PNG24:".sprintf("./images/training_gen/L_Hx%02dHy%02d_Vx%02dVy%02d",$x,$y+$size,$x,$y).".png");


		$img = $imgTmp1->clone();
		#horizontal for T
		$img->Draw(primitive=>'line',stroke=>'black', strokewidth=>'1', fill=>'none', points=>($x).",".($y). " ".($x+$size).",".($y));
		#vertical for T		
		$img->Draw(primitive=>'line',stroke=>'black', strokewidth=>'1', fill=>'none', points=>($x+$size/2).",".($y). " ".($x+$size/2).",".($y+$size));

		$img->Write("PNG24:".sprintf("./images/training_gen/T_Hx%02dHy%02d_Vx%02dVy%02d",$x,$y,$x+$size/2,$y).".png");


		
	}
}
