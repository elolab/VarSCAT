def Cfunc(df):
	C_bed = open("Complex_variants_bed.bed","a")
	for i in range(0,df.shape[0]):
		if df.loc[i,"UPS_LEFT"] != "":
			index.append(i)
			item = str(df.loc[i,"Chromosome"])+"	"+str(df.loc[i,"UPS_L"])+"	"+str(df.loc[i,"UPS_R"])+"\n"
			C_bed.write(item)
	
	C_bed.close()
	
def Rfunc(df):
	R_bed = open("satellite_DNA_region_variants.bed","a")
	for i in range(0,df.shape[0]):
		item = str(df.loc[i,"Chromosome"])+"	"+str(df.loc[i,"Start"])+"	"+str(df.loc[i,"End"])+"\n"
		R_bed.write(item)
	
	R_bed.close()
