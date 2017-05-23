import pandas as pd
#file=pd.read_csv('',sep='\t')
#f7=file[file['Tail Chromosome']=='chr7']
#factb=f7[f7['Tail Start'].between(5525148,5532601)]
#factb.to_csv('')
factb=pd.read_csv('')

#filter by AU content?
def aucontent(st):
	l=len(st)
	if st.count('A')>l*0.5 or st.count('T')>l*0.5:
		return True
	else:
		return False
def austring(st):
	l=len(st)
	if 'AAA' in st or 'TTT' in st:
		return True
	else:
		return False
filtered0=f0[f0.apply(lambda x: aucontent(x['Tail']), axis=1)]
filtered2=f2[f2.apply(lambda x: aucontent(x['Tail']), axis=1)]
#make bedgraph file
he='''track type=bedGraph name="Tail Counts" graphType=bar\n'''
with open('2-2f.bedgraph','w') as f:
	f.write(he)
newf=pd.DataFrame()
newf['chr']=factb['Tail Chromosome']
newf['start']=factb['Tail Start']
newf['stop']=newf['start']+1
newf['count']=1
newf=newf.groupby(['chr','start','stop']).sum()
newf.to_csv('2-2f.bedgraph',sep='\t',header=False,mode='a')
