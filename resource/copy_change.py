from .repeatAlign import diff_letters

def p_copyN(v_ref,v_alt,v_start,v_end,r_start,r_end,r_motif,r_period):	
	cp_change = "#"
	if len(v_ref) > len(v_alt):
		if v_start <= r_start and (v_start+len(v_ref)-len(v_alt)) >=r_end:
			cp_change = "-"+str(r_period)
		elif v_start > r_start and (v_start+len(v_ref)-len(v_alt)) < r_end:
			v_pattern = v_ref
			if len(v_pattern)>len(r_motif):
				cp=0
				while len(v_pattern)>len(r_motif):
					if v_pattern[:len(r_motif)]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[len(r_motif):]
					else:
						if cp==0:
							v_pattern=v_pattern[1:]
						else:
							break
			elif len(v_pattern)==len(r_motif):
				cp=0
				if v_pattern==r_motif:
					cp=cp+1
				else:
					n=0
					while n<len(r_motif):
						r_motif=r_motif[1:]+r_motif[0]
						if v_pattern==r_motif:
							cp=cp+1
							break
						n=n+1
			else:
				cp=0
				
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)
		elif v_start >= r_end:
			cp_change = "=0"
		elif v_end<r_start:
			cp_change = "=0"
		elif v_start == r_start:
			v_pattern = v_ref
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[-len(r_motif):]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[:len(v_pattern)-len(r_motif)]
					else:
						break
			else:
				cp=0	
			
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)
		elif v_start < r_start and (v_start+len(v_ref)-len(v_alt)) >= r_start and v_start+len(v_ref)-len(v_alt) < r_end:
			v_pattern = v_ref[r_start-v_end-2:len(v_ref)]
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[-len(r_motif):]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[:len(v_pattern)-len(r_motif)]
					else:
						if cp==0:
							v_pattern=v_pattern[:-1]
						else:
							break
			elif len(v_pattern)<len(r_motif):
				cp=0
		
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)
		elif v_start > r_start and v_start <= r_end and (v_start+len(v_ref)-len(v_alt)) > r_end:
			v_pattern = v_ref[0:r_end-v_start+2]
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[:len(r_motif)]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[len(r_motif):]
					else:
						if cp==0:
							v_pattern=v_pattern[1:]
						else:
							break
			elif len(v_pattern)<len(r_motif):
				cp=0		
		
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)
		else:	
			v_pattern=v_ref[len(v_alt):]
			cp=d_change_inner(v_pattern,r_motif,0,0)
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)	
	elif len(v_ref) < len(v_alt):
		if v_start >= r_start and v_end <= r_end:
			v_pattern=v_alt[len(v_ref):]
			if len(v_pattern)>len(r_motif):
				cp=0
				while len(v_pattern)>len(r_motif):
					if v_pattern[:len(r_motif)]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[len(r_motif):]
					else:
						break
			elif len(v_pattern)<=len(r_motif):
				cp=0
				
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "+"+str(cp)
		elif v_start >= r_end:
			cp_change = "=0"
		elif v_start == r_start -1:
			v_pattern = v_alt
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[-len(r_motif):]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[:len(v_pattern)-len(r_motif)]
					else:
						break
			else:
				cp=0	
			
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "+"+str(cp)
		elif v_end == r_start-1 :
			v_pattern = v_alt
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[-len(r_motif):]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[:len(v_pattern)-len(r_motif)]
					else:
						break
			else:
				cp=0
				
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "+"+str(cp)		
		else:
			v_pattern=v_alt[len(v_ref):]
			cp=d_change_inner(v_pattern,r_motif,0,0)
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "+"+str(cp)	
	elif len(v_ref) == len(v_alt):
		if v_start >= r_start and v_end <= r_end:
			v_pattern = v_ref
			if len(v_pattern)>len(r_motif):
				cp=0
				while len(v_pattern)>len(r_motif):
					if v_pattern[:len(r_motif)]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[len(r_motif):]
					else:
						if cp==0:
							v_pattern=v_pattern[1:]
						else:
							break
			elif len(v_pattern)==len(r_motif):
				cp=0
				if v_pattern==r_motif:
					cp=cp+1
				else:
					n=0
					while n<len(r_motif):
						r_motif=r_motif[1:]+r_motif[0]
						if v_pattern==r_motif:
							cp=cp+1
							break
						n=n+1
			else:
				cp=1
				
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)
		elif v_end < r_start:
			v_pattern = v_alt
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[-len(r_motif):]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[:len(v_pattern)-len(r_motif)]
					else:
						break
			elif len(v_pattern)<len(r_motif):
				cp=0
				
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "+"+str(cp)
		elif v_start > r_end:
			v_pattern = v_alt
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[:len(r_motif)]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[len(r_motif):]
					else:
						break
			elif len(v_pattern)<len(r_motif):
				cp=0		
			
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "+"+str(cp)
		elif v_start < r_start and v_end >= r_start and v_end < r_end:
			v_pattern = v_ref[r_start-v_end-2:len(v_ref)]
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[-len(r_motif):]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[:len(v_pattern)-len(r_motif)]
					else:
						break
			elif len(v_pattern)<len(r_motif):
				cp=0
		
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)
		
		elif v_start > r_start and v_start <= r_end and v_end > r_end:
			v_pattern = v_ref[0:r_end-v_start+2]
			if len(v_pattern)>=len(r_motif):
				cp=0
				while len(v_pattern)>=len(r_motif):
					if v_pattern[:len(r_motif)]==r_motif:
						cp=cp+1
						v_pattern=v_pattern[len(r_motif):]
					else:
						break
			elif len(v_pattern)<len(r_motif):
				cp=0		
		
			if cp==0:
				cp_change = "=0"
			else:
				cp_change = "-"+str(cp)
		elif v_start < r_start and v_end > r_end:
			cp_change = "-"+str(r_period)
		else:
			cp_change = "=0"
	return cp_change

def d_change_inner(v_pattern1,r_motif1,r_tolerate1,g_tolerate1):
	cp1=0
	gap_n=0
	if len(v_pattern1)==len(r_motif1):
		if v_pattern1==r_motif1:
			cp1=1
		else:
			n=0
			while n<len(r_motif1):
				r_motif1=r_motif1[1:]+r_motif1[0]
				if v_pattern1==r_motif1:
					cp1=1
					break
				n=n+1
				
	elif len(v_pattern1)>len(r_motif1):
		while len(v_pattern1)>=len(r_motif1):
			pattern_check=v_pattern1[0:len(r_motif1)]		
			common_base = diff_letters(pattern_check,r_motif1)
			if float(common_base/len(r_motif1))>=r_tolerate1:
				cp1=cp1+1
				v_pattern1=v_pattern1[len(r_motif1):]
				gao_n=0
			else:
				gap_n=gap_n+1
				if g_tolerate1 == -1:
					g_tolerate1 = len(r_motif1)
				if gap_n<=g_tolerate1:
					v_pattern1=v_pattern1[1:]
				else:
					break
	return cp1

def d_copyN(v_ref,v_alt,v_start,v_end,r_start,r_end,r_motif,r_period,r_tolerate,g_tolerate):	
	if len(v_ref)>len(v_alt):
		v_pattern=v_ref[len(v_alt):]
		if v_start<=r_start and v_start+len(v_pattern)-1>=r_end:
			cp=r_period
		else:
			cp=d_change_inner(v_pattern,r_motif,r_tolerate,g_tolerate)
	elif len(v_ref)<len(v_alt):
		v_pattern=v_alt[len(v_ref):]
		cp=d_change_inner(v_pattern,r_motif,r_tolerate,g_tolerate)
	elif len(v_ref)==len(v_alt):
		v_pattern=v_ref
		if v_start<=r_start and v_end>=r_end:
			cp=r_period
		else:
			cp=d_change_inner(v_pattern,r_motif,r_tolerate,g_tolerate)
	
	cp_change = "#"
	if cp == 0:
		cp_change = "=0"
	else:
		if len(v_ref)>len(v_alt):
			cp_change = "-"+str(cp)
		elif len(v_ref)<len(v_alt):
			cp_change = "+"+str(cp)
		elif len(v_ref)==len(v_alt):
			cp_change = "-"+str(cp)
	return cp_change

def check_copy(v_ref0,v_alt0,v_start0,v_end0,r_start0,r_end0,r_motif0,r_period0,match_score,r_tolerate,g_tolerate):	
	if match_score==1:
		copy_number_change=p_copyN(v_ref0,v_alt0,v_start0,v_end0,r_start0,r_end0,r_motif0,r_period0)
	else:
		copy_number_change=d_copyN(v_ref0,v_alt0,v_start0,v_end0,r_start0,r_end0,r_motif0,r_period0,r_tolerate,g_tolerate)
	
	return copy_number_change




