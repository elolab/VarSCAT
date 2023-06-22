import pybedtools
import os

def annotate_region(bedin_input,df_input):
	location_list=[]
	info_list=[]
	for i in range(0,len(df_input)):
		variant_r = [(df_input[i][0], df_input[i][4], df_input[i][6]),]
		v_bed = pybedtools.BedTool(variant_r)
		ann = v_bed.intersect(bedin_input,wa=True,wb=True)
		location_string = ""
		info_string = ""
		for j in range(0,len(ann)):
			location_string=location_string+ann[j][3]+":"+ann[j][4]+"-"+ann[j][5]
			for m in range(6,len(ann[j].fields)):
				if m!=len(ann[j].fields)-1:
					info_string=info_string+ann[j][m]+"|"
				else:
					info_string=info_string+ann[j][m]
			if j<(len(ann)-1) and len(ann)!=1:
				location_string=location_string+";"
				info_string=info_string+";"
		location_list.append(location_string)
		info_list.append(info_string)
	
	return location_list,info_list

def custom_anno(df,loc,bed_file,anno_file):
	if not os.path.isabs(anno_file):
		anno_file=os.path.normpath(os.path.join(os.getcwd(),anno_file))
	a = pybedtools.example_bedtool(anno_file)
	
	if loc!="":
		chromosome=loc.split(":")[0]
		coordinates=loc.split(":")[1]
		start=int(coordinates.split("-")[0])
		end=int(coordinates.split("-")[1])
		location = [(chromosome, start, end),]
		b = pybedtools.BedTool(location)
		bedin = a.intersect(b,wa=True)
		anno_results=annotate_region(bedin,df)
	elif bed_file!="":
		if not os.path.isabs(bed_file):
			bed_file=os.path.normpath(os.path.join(os.getcwd(),bed_file))
		b = pybedtools.example_bedtool(bed_file)
		bedin = a.intersect(b,wa=True)
		anno_results=annotate_region(bedin,df)

	return anno_results

