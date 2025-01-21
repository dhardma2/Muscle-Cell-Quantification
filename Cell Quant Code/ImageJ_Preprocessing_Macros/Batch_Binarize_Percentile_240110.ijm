//Macro for preprocessing and binarising static images stained for nuclei and
// myotubes with the option to include myogenic specific stainings and other 
//stainings (such as BTX)
  
Dir = getDirectory("Choose a Directory ");

File_struct = "Day4 Untrained CTR-01-Stitching-";
file_type=".czi";
Dialog.create("New Image");
Dialog.addString("File name structure:", File_struct);
Dialog.addString("File Type:", file_type);
Dialog.addString("Experiment ID:", "Untrained");
Dialog.addNumber("Number of images:", 2);
Dialog.addChoice("All nuclei channel:", newArray("1", "2", "3", "4","N/A"));
Dialog.addNumber("All nuclei  threshold:", 90);
Dialog.addChoice("Myonuclei channel:", newArray("1", "2", "3", "4","N/A"));
Dialog.addNumber("Myonuclei  threshold:", 95);
Dialog.addChoice("Mytube channel:", newArray("1", "2", "3", "4","N/A"));
Dialog.addNumber("Myotube  threshold:", 75);
Dialog.addChoice("Extra staining channel:", newArray("1", "2", "3", "4","N/A"));
Dialog.addString("Extra staining suffix:", "BTX");
Dialog.addNumber("Extra staining threshold:", 75); 
Dialog.addChoice("z-projection:", newArray("Yes", "No"));
Dialog.addChoice("Batch mode on:", newArray("True", "False"));
Dialog.show();
File_struct = Dialog.getString();
file_type = Dialog.getString();
Expt= Dialog.getString();
no_imgs = Dialog.getNumber();
DAPI_ch=Dialog.getChoice();
DAPI_thresh= Dialog.getNumber();
MTDNA_ch= Dialog.getChoice();
MTDNA_thresh= Dialog.getNumber();
MT_ch= Dialog.getChoice();
MT_thresh= Dialog.getNumber();
Other_ch= Dialog.getChoice();
Other_suff= Dialog.getString();
Other_thresh= Dialog.getNumber();
z_proj= Dialog.getChoice();
Batch_mode= Dialog.getChoice();

//pixel size used for rolling ball background subtraction
SubtractionPxSize="70";

for (i = 1; i < no_imgs; i++) {
// select file
if(i < 10){
	filename=File_struct+"0"+i;
}
else filename=File_struct+i;
//apply batch mode if required
if (Batch_mode == "True"){
	setBatchMode(true);
}
//import images
Importer="open=["+Dir+filename+file_type+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2";
run("Bio-Formats Importer",Importer );
if (z_proj == "Yes"){
//z-stack based on average intensity
	run("Z Project...", "projection=[Average Intensity]");
//background subtraction
	run("Subtract Background...", "rolling="+SubtractionPxSize+" stack");
}
rename("Pos"+i);
Pos=getTitle;
run("Split Channels");

function Percentile_Threshold(percentage){
nBins = 256; 
resetMinAndMax(); 
getHistogram(values, counts, nBins); 
// find culmulative sum 
nPixels = 0; 
for (i = 0; i<counts.length; i++) 
  nPixels += counts[i]; 
nBelowThreshold = nPixels * percentage / 100; 
sum = 0; 
for (i = 0; i<counts.length; i++) { 
  sum = sum + counts[i]; 
  if (sum >= nBelowThreshold) { 
    setThreshold(values[0], values[i]); 
    print(values[0]+"-"+values[i]+": "+sum/nPixels*100+"%"); 
    i = 99999999;//break 
  } 
} 
}

function Make_Binary(im_title,percentage){
saveAs("Tiff", im_title+".tif");
Percentile_Threshold(percentage);
run("Convert to Mask");
run("Invert");
saveAs("Tiff", im_title+"BW.tif");             
close();

}
selectWindow("C"+DAPI_ch+"-"+Pos);
title=Dir+Expt+"/"+Expt+"-"+Pos+"-DNA"; 
Make_Binary(title,DAPI_thresh);
selectWindow("C"+MTDNA_ch+"-"+Pos);
if (MTDNA_ch != "N/A"){
	title=Dir+Expt+"/"+Expt+"-"+Pos+"-MTDNA";
	Make_Binary(title,MTDNA_thresh);
}
selectWindow("C"+MT_ch+"-"+Pos);
title=Dir+Expt+"/"+Expt+"-"+Pos+"-Actin";
Make_Binary(title,MT_thresh);
if (Other_ch != "N/A"){
	selectWindow("C"+Other_ch+"-"+Pos);	
	title=Dir+Expt+"/"+Expt+"-"+Pos+"-"+Other_suff;
	Make_Binary(title,Other_thresh);
}
};
