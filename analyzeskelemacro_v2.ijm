//output_dir = "Y:/fast/_SR/si/cwolf/Analysis Clips/output.COMP/";
output_dir = "Z:/_SR/si/cwolf/Analysis Clips/output.COMP/";
name=getTitle();


run("Split Channels");
selectImage(name+" (blue)");
close();
selectImage(name+" (green)");
close();

selectImage(name+" (red)");

name=getTitle();
name = substring(name,0,lengthOf(name)-4);
rename(name);
//if(selectionType != -1) run("Clear Outside");
setOption("BlackBackground", false);
run("Convert to Mask");
run("Fill Holes");
run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
run("Despeckle");
run("Despeckle");
run("Despeckle");
run("Skeletonize");
run("Analyze Skeleton (2D/3D)", "prune=[lowest intensity voxel] prune_0 calculate show display original_image=["+name+"]");
saveAs("PNG", output_dir+name+"-labeled-skeletons.png");
close();

saveAs("Results", output_dir+name+"_RESULTS.csv");
close("Longest shortest paths");
close("Results");
close(name);
close("Branch information");