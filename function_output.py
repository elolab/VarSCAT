def Cfunc(df,output,sample):
	name=output+"_ambiguity_variants.bed"
	C_bed = open(name,"a")
	item="#Chromosome"+"	"+"Ambiguity_left"+"	"+"Ambiguity_right"+"	"+"Position"+"	"+"REF"+"	"+"ALT"
	for j in range(0,len(sample)):
		item=item+"	"+sample[j]
	item=item+"\n"
	C_bed.write(item)
	for i in range(0,df.shape[0]):
		if df.loc[i,"UPS_LEFT"] != "":
			item = str(df.loc[i,"Chromosome"])+"	"+str(df.loc[i,"UPS_L"])+"	"+str(df.loc[i,"UPS_R"])+"	"+str(df.loc[i,"Position"])+"	"+str(df.loc[i,"REF"])+"	"+str(df.loc[i,"ALT"])
			for j in range(0,len(sample)):
				item=item+"	"+str(df.loc[i,sample[j]])
			item=item+"\n"
			C_bed.write(item)
	
	C_bed.close()
	
def Rfunc(df,output):
	name=output+"_TR_region_variants.bed"
	R_bed = open(name,"a")
	item="#Chromosome"+"	"+"Repeat_Start"+"	"+"Repeat_End"+"	"+"Motifs"+"	"+"Period"+"	"+"Size"+"	"+"Repeat_GC%"+"	"+"Repeat_Score"+"	"+"Alignment_Score"+"	"+"Gap_Score"+"\n"
	R_bed.write(item)
	for i in range(0,df.shape[0]):
		if df.loc[i,"Motifs"] != "":
			if str(df.loc[i,"Start"]).find(",") == -1:
				item = str(df.loc[i,"Chromosome"])+"	"+str(df.loc[i,"Start"])+"	"+str(df.loc[i,"End"])+"	"+str(df.loc[i,"Motifs"])+"	"+str(df.loc[i,"Period"])+"	"+str(df.loc[i,"Size"])+"	"+str(df.loc[i,"Repeat_GC%"])+"	"+str(df.loc[i,"Repeat_Score"])+"	"+str(df.loc[i,"Alignment_Score"])+"	"+str(df.loc[i,"Gap_Score"])+"\n"
				R_bed.write(item)
			else:
				start=str(df.loc[i,"Start"]).split(",")
				end=str(df.loc[i,"End"]).split(",")
				motifs=str(df.loc[i,"Motifs"]).split(",")
				period=str(df.loc[i,"Period"]).split(",")
				size=str(df.loc[i,"Size"]).split(",")
				GC=str(df.loc[i,"Repeat_GC%"]).split(",")
				RS=str(df.loc[i,"Repeat_Score"]).split(",")
				AS=str(df.loc[i,"Alignment_Score"]).split(",")
				GS=str(df.loc[i,"Gap_Score"]).split(",")
				for i in range(0,len(start)):
					item = str(df.loc[i,"Chromosome"])+"	"+start[i]+"	"+end[i]+"	"+motifs[i]+"	"+period[i]+"	"+size[i]+"	"+GC[i]+"	"+RS[i]+"	"+AS[i]+"	"+GS[i]+"\n"
					R_bed.write(item)
	
	R_bed.close()

