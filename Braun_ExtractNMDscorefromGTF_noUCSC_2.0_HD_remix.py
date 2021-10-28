import sys
import re
import argparse
import pandas as pd
import numpy as np
import numbers
from gtfparse import read_gtf
import subprocess
import gzip
from Bio import SeqIO
from Bio.Seq import Seq

##Abre el archivo Braun.csv
filename = '/home/vant/Braun/Braun.csv'

## Carga archivos de referencia (canon, gtf, hgid)
# canon = pd.read_csv('/home/vant/NMDetective/canonical.txt', sep='\t')
# canon_init = canon.copy() # keep input canon here for using it in getPTCpositionFromFS function
# canon = canon[['transcript']]
# # Otra forma de filtrar es hacer un inner join con la tabla de canonicos. Nos quedamos con esta opción porque canon sale de UCSC y tambien lo usamos para sacar los id UCSC del GTF.
# canon.columns = ['Feature']
# canon_aux = canon.copy()
# canon_aux['Feature'] = canon_aux['Feature'].str.split('.').str[0]
gtf = read_gtf("/home/vant/NMDetective/CDS_hg19_NMDetectiveB_Lindeboom_et_al.v2.gtf")
gtf_all = read_gtf("/home/vant/NMDetective/hg19_NMDetectiveB_Lindeboom_et_al.v2.gtf")
gtf_all = gtf_all.loc[(gtf_all['feature'] == 'exon') | (gtf_all['feature'] == 'CDS')]

# Abrir las tablas necesarias, la del ENST y la de conónicos
# kgid = pd.read_csv('/home/vant/NMDetective/kgAlias.txt', sep='\t', header=None)
# # Poner nombre a las columnas de kgAlias.txt
# kgid.columns = ['Feature', 'gene_id']

# Aqui GTF solo contiene los canonicos
# gtf = pd.merge(gtf, kgid, on='gene_id')  # Para unir el KgID al ENST
# gtf = pd.merge(gtf, canon, on='Feature')  # Para selecionar el trancrito canónico
gtf['Feature'] = gtf['exon_id'].str.split('.').str[0]
# Genoma de referencia
refgen = '/home/vant/Braun/ucsc.hg19.fasta'

def nmdxtract(ID, PTCposition, mega_gtf, consequence):
	##Sacamos los cromosomas y las posiciones
	CHR = ID.split('_')[0]
	if 'Nonsense' in consequence:
		POS = int(ID.split('_')[1])
		mega_gtf = mega_gtf.loc[(mega_gtf['seqname'] == CHR) & (mega_gtf['start'] <= POS) & (mega_gtf['end'] >= POS)]
		if mega_gtf.empty == False:
			start = list(mega_gtf['start'].values)[0]
			end = list(mega_gtf['end'].values)[0]
			nmd = mega_gtf['nmd_score'].head(1).values[0].split(',')
			strand = str(mega_gtf['strand'])
			r = np.arange(start, end + 1, 1)
			r = r.tolist()
			d = dict(zip(r, nmd))
			return d[POS]
		else:
			return '.'
	else:
		if (PTCposition != '.') & (PTCposition != 'utr') & (PTCposition != 'out'):
			POS = int(PTCposition)
			mega_gtf = mega_gtf.loc[(mega_gtf['seqname'] == CHR) & (mega_gtf['start'] <= POS) & (mega_gtf['end'] >= POS)]
			if mega_gtf.empty == False:
        	##Extraemos del GTF las filas que coinciden (tendría que salir 1, mirar esto)
        	#mega_gtf.loc[(mega_gtf['seqname'] == CHR) & (mega_gtf['start'] <= POS) & (mega_gtf['end'] >= POS)]
        	##Buscar el primer número del intervalo
				start = list(mega_gtf['start'].values)[0]
				##Buscar el segundo número del intervalo
				end = list(mega_gtf['end'].values)[0]
				##Extrae los valores de NMD del intervalo y lo transforma a lista separando por comas
				nmd = mega_gtf['nmd_score'].head(1).values[0].split(',')
				# Extrae la strand para ver si es forward o reverse
				strand = str(mega_gtf['strand'])
				##Crea el intervalo, se pone +1 para que coincida con len(nmd)
				r = np.arange(start, end + 1, 1)
				r = r.tolist()
				##Crea un diccionario que une posición en el intervalo con score de NMD
				d = dict(zip(r, nmd))
				return d[POS]
			else:
				return '.'
		elif PTCposition  == '.':
			return '.'
		elif PTCposition == 'utr':
			return '0.00'
		else:
			return '0.666'

def getPTCpositionFromFS(gtf, ID, VCF_init, refgen, consequence, ID_specific):
	if 'Nonsense' in consequence:
		return '.'
	else:
		ID = ID.split('_')
		CHR = ID[0]
		POS = int(ID[1])
		REF = ID[2]
		ALT = ID[3]
		ID_aux = '_'.join([CHR,str(POS)])

		if gtf.loc[(gtf['seqname'] == CHR) & (gtf['start'] <= POS) & (gtf['end'] >= POS)].empty == False:
			feature = gtf['Feature'].loc[(gtf['seqname'] == CHR) & (gtf['start'] <= POS) & (gtf['end'] >= POS)].values[0]
			exons = gtf.loc[gtf['exon_id'].str.contains(feature)]
			strand = exons['strand'].values[0]
			ucsc_id = exons['gene_id'].values[0]
			exons['exon_number'] = exons['exon_number'].astype(int)
			if strand == '+': #Para añadir zona UTR
				if gtf_all['exon_number'].values[-1] == exons['exon_number'].values[-1]:
					gtf_aux = gtf_all.loc[gtf_all['gene_id'] == ucsc_id]
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['start'] == (exons['start'].values[-1]))]
					row_utr['start'] = (exons['end'].values[-1]) + 1
					row_utr['exon_number'] = row_utr['exon_number'].astype(int)
					row_utr['exon_number'] = row_utr['exon_number'] + 1
					exons = pd.concat([exons, row_utr], ignore_index = True, axis = 0)
				else:
					gtf_aux = gtf_all.loc[gtf_all['gene_id'] == ucsc_id]
					gtf_aux['exon_number'] = gtf_aux['exon_number'].astype(int)
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['exon_number'] > (exons['exon_number'].values[-1]))]
					exons = pd.concat([exons, row_utr], ignore_index = True, axis = 0)
			else:
				if gtf_all['exon_number'].values[0] == exons['exon_number'].values[0]:
					gtf_aux = gtf_all.loc[gtf_all['gene_id'] == ucsc_id]
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['end'] == (exons['end'].values[0]))]
					row_utr['end'] = (exons['start'].values[0]) - 1
					row_utr['exon_number'] = row_utr['exon_number'].astype(int)
					row_utr['exon_number'] = row_utr['exon_number'] - 1
					exons = pd.concat([row_utr, exons], ignore_index = True, axis = 0)
				else:
					gtf_aux = gtf_all.loc[gtf_all['gene_id'] == ucsc_id]
					gtf_aux['exon_number'] = gtf_aux['exon_number'].astype(int)
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['exon_number'] < (exons['exon_number'].values[0]))]
					exons = pd.concat([row_utr, exons], ignore_index = True, axis = 0)
			#if strand == '-': #probar quitandolo para ver si funciona
				#exons['exon_number'] = list(exons['exon_number'])[::-1]
			mutExon = exons.loc[(exons['seqname'] == CHR) & (exons['start'] <= POS) & (exons['end'] >= POS)]
			exons['ID'] = exons['seqname'] + ':' + exons['start'].astype(str) + '-' + exons['end'].astype(str)
			filein = ID_aux + '.exons.csv'
			exons['ID'].to_csv(filein, sep='\t', index = None, header = None)
			#Get nucleotide sequences from all exons
			fileout = ID_aux + '.outExons.fa'
			command = ' '.join(['samtools faidx', refgen, '--region-file', filein, '>', fileout])
			subprocess.run(command, shell=True)
			exonsSeqsdict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fileout, "fasta")}

			# Get a vcf containing only the variant of interest
			filein = ID_aux + '.mutation.txt'
			f = open(filein, 'w')
			f.write('\t'.join(ID))
			f.close()
			fileout = ID_aux + '.mut.vcf'
			command = ' '.join(['zgrep  ^#', VCF_init, '>', fileout, '; zgrep ', ID_specific, VCF_init, '>>', fileout, '; bgzip -c', fileout, '>', fileout +'.gz','; tabix -p vcf', fileout +'.gz'])
			subprocess.run(command, shell=True)
			fileout = ID_aux + '.mut.vcf.gz'

			# Get mutated sequence from the specific exon
			filein = ID_aux + '.mutExon.csv'
			(mutExon['seqname'] + ':' + mutExon['start'].astype(str) + '-' + mutExon['end'].astype(str)).to_csv(filein,sep='\t',index = None, header = None)
			filein2 = fileout
			fileout = ID_aux + '.outmutExon.fa'
			command = ' '.join(['samtools faidx', refgen, '--region-file', filein, '| bcftools consensus', filein2, '-o', fileout]) #-s TUMOR
			subprocess.run(command, shell=True)
			mutExonsSeqsdict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fileout, "fasta")}
			exonsSeqsdict.update(mutExonsSeqsdict)
			exons['sequence'] = exons['ID'].map(exonsSeqsdict)
			exons['exon_number'] = exons['exon_number'].astype(int)
			# offset_fs = (len(ALT.replace('-', '')) - len(REF.replace('-', ''))) % 3
			# offset_utr = len(exouc002lyj.2ns['sequence'].loc[exons['feature'] == 'exon']) % 3
			# offset = offset_fs + offset_utr
			exons = exons.sort_values('exon_number', ascending=True)


			if strand == '+':
				exons['ExonLen'] = exons['sequence'].apply(lambda x: len(str(x))) #Sacar lenght usando la sequence porque si usas exon start - end vas perdiendo 1 en cada exón
				exons['translen'] = exons['ExonLen'].cumsum()
				exons['translenplusone'] = exons['translen'] + 1
				exons['startTranscript'] = exons['translenplusone'].shift(periods=1)
				exons['startTranscript']  = exons['startTranscript'].fillna(0).astype(int)
				startCDS = exons['start'].values[0]
				offset = len(''.join(list(exons['sequence']))) % 3
				sequence = Seq(''.join(list(exons['sequence']))[:-offset])
				#mutcodon = int((np.floor((POS - (startCDS)) /3)) + 1)
				#mutcodonposition = (POS - startCDS) % 3

				sequence_PTC = str(sequence.translate(to_stop=True))

				PTCposinTrans = len(sequence_PTC)*3
				if PTCposinTrans >= exons['translen'].max():
					return 'out'
				else:
					if 'Del' in consequence:
						#PTCposinTransCorrection = PTCposinTrans + (len(REF)-1)
						exons_aux = exons.loc[(exons['startTranscript'] <= PTCposinTrans) & (exons['translenplusone'] > PTCposinTrans)]
						if exons_aux.loc[exons_aux['feature'] == 'exon'].empty == False:
							return 'utr'
						else:
							Distance_to_start_exon = PTCposinTrans - exons_aux['startTranscript'].values[0]
							PTCposition = exons_aux['start'].values[0] + Distance_to_start_exon + (len(REF)-1)
							return PTCposition
					else:
						#PTCposinTransCorrection = PTCposinTrans - (len(ALT)-1)
						exons_aux = exons.loc[(exons['startTranscript'] <= PTCposinTrans) & (exons['translenplusone'] > PTCposinTrans)]
						if exons_aux.loc[exons_aux['feature'] == 'exon'].empty == False:
							return 'utr'
						else:
							Distance_to_start_exon = PTCposinTrans - exons_aux['startTranscript'].values[0]
							PTCposition = exons_aux['start'].values[0] + Distance_to_start_exon - (len(ALT)-1)
							return PTCposition

			else:
				exons['ExonLen'] =  exons['sequence'].apply(lambda x: len(str(x)))
				exons['translen'] = exons.loc[::-1, 'ExonLen'].cumsum()[::-1]
				exons['translenplusone'] = exons['translen'] + 1
				exons['startTranscript'] = exons['translenplusone'].shift(periods=-1)
				exons['startTranscript']  = exons['startTranscript'].fillna(0).astype(int)
				startCDS = exons.sort_values('exon_number', ascending=False)['end'].values[0]
				offset = len(''.join(list(exons['sequence']))) % 3
				sequence = Seq(''.join(list(exons.sort_values('exon_number', ascending=True)['sequence']))[offset:]).reverse_complement()
				#mutcodon = int((np.floor((startCDS - POS) / 3)) + 1)
				#mutcodonposition = ((startCDS) - POS) % 3
				sequence_PTC = str(sequence.translate(to_stop=True))
				PTCposinTrans = len(sequence_PTC)*3
				if PTCposinTrans >= exons['translen'].max():
					return 'out'
				else:
					if 'Del' in consequence:
						#PTCposinTransCorrection = PTCposinTrans + (len(REF)-1)
						exons_aux = exons.loc[(exons['startTranscript'] <= PTCposinTrans) & (exons['translenplusone'] > PTCposinTrans)]
						if exons_aux.loc[exons_aux['feature'] == 'exon'].empty == False:
							return 'utr'
						else:
							Distance_to_start_exon = PTCposinTrans - exons_aux['startTranscript'].values[0]
							PTCposition = exons_aux['end'].values[0] - Distance_to_start_exon + (len(REF)-1)
							return PTCposition
					else:
						#PTCposinTransCorrection = PTCposinTrans - (len(ALT)-1)
						exons_aux = exons.loc[(exons['startTranscript'] <= PTCposinTrans) & (exons['translenplusone'] > PTCposinTrans)]
						if exons_aux.loc[exons_aux['feature'] == 'exon'].empty == False:
							return 'utr'
						else:
							Distance_to_start_exon = PTCposinTrans - exons_aux['startTranscript'].values[0]
							PTCposition = exons_aux['end'].values[0] - Distance_to_start_exon - (len(ALT)-1)
							return PTCposition
		else:
			return '.'

# Read Braun's mutation file
# OUTPUT list stored into pandas df


df = pd.read_csv(filename, sep='\t')
REF = 'Tumor_Seq_Allele1'
ALT = 'Tumor_Seq_Allele2'
chr = 'Chromosome'
position = 'Start_position'

# Filtras para quedarte con las PTC en los trasncript canonicos del vcf
df = df.loc[df['Variant_Classification'].str.contains('Nonsense|Frame_Shift')]
#idx = (df['Variant_Classification'] == 'Frame_Shift_Ins')
#df_aux = df.copy()
#df_aux.loc[idx,['Start_position','End_position']] = df_aux.loc[idx,['End_position','Start_position']].values
#df_aux['VEP_allele'] = df_aux['Tumor_Seq_Allele1'] + '/' + df_aux['Tumor_Seq_Allele2']
#df_aux['ID_specific'] = df_aux['Chromosome'].astype(str) + '_' + df_aux['Start_position'].astype(str)  + '_' + df_aux['VEP_allele'] + '_' + df_aux['Tumor_Sample_Barcode']
#df['ID_specific'] = df['ID'].map(dict(zip(list(df_aux['ID']),list(df_aux['ID_specific']))))

df['VEP_allele'] = df[REF] + '/' + df[ALT]
df['ID_specific'] = df[chr].astype(str) + '_' + df[position].astype(str) + '_' + df['VEP_allele'] + '_' + df['Tumor_Sample_Barcode']
df[chr] = 'chr' + df[chr].astype(str)

########## MODIFICATION OF INDEL VARIANTS IN AGREEMENT WITH VCF_init
# Insertions, adding the reference
def addREF2Insertions(CHR, POS, refgen):
    region = CHR + ':' + str(POS) + '-' + str(POS)
    fileout = CHR + '_' + str(POS) + '_' + str(POS) + '.fa'
    command = ' '.join(['samtools faidx', refgen, region, '>', fileout])
    subprocess.run(command, shell=True)
    Seqsdict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fileout, "fasta")}
    command = ' '.join(['rm', fileout])
    subprocess.run(command, shell=True)
    return Seqsdict[region]

df['aux'] = df[REF]
df[REF] = df.apply(lambda x: addREF2Insertions(x[chr], x[position], refgen).upper() if x['aux'] == '-' else x[REF], axis = 1)
df[ALT] = df.apply(lambda x:  x[REF] + x[ALT] if x['aux'] == '-' else x[ALT], axis=1)

# Deletions, adding the previous position as reference
def addREF2deletions(CHR, POS, refgen):
    region = CHR + ':' + str(POS) + '-' + str(POS)
    fileout = CHR + '_' + str(POS) + '_' + str(POS) + '.fa'
    command = ' '.join(['samtools faidx', refgen, region, '>', fileout])
    subprocess.run(command, shell=True)
    Seqsdict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fileout, "fasta")}
    command = ' '.join(['rm', fileout])
    subprocess.run(command, shell=True)
    return Seqsdict[region]

df['aux'] = df[ALT]
# Position -1 for deletions
df[position] = df.apply(lambda x: x[position] - 1 if x['aux'] == '-' else x[position], axis=1)
df[ALT] = df.apply(lambda x: addREF2deletions(x[chr], x[position], refgen).upper() if x['aux'] == '-' else x[ALT], axis = 1)
df[REF] = df.apply(lambda x: x[ALT] + x[REF] if x['aux'] == '-' else x[REF], axis = 1)

df.drop(columns='aux',inplace = True)
####################################################################
####################################################################
df['ID'] = df[chr].astype(str) + '_' + df[position].astype(str) + '_' + df[REF] + '_' + df[ALT]

#VCF_init = '/mnt/64716603-5b56-4f9a-b195-c11560647a3a/data/TFM_NMD/Braun/sorted.chr.Annot_input_VEP_Braun.vcf.gz'
VCF_init = '/home/vant/Braun/Braun.vcf.gz'
df['PTCposition'] = df.apply(lambda x: getPTCpositionFromFS(gtf, x['ID'], VCF_init, refgen, x['Variant_Classification'], x['ID_specific']), axis=1)

df = df.reset_index()
df.drop('index',1, inplace = True)


# NMD = []
# to_5 = []
# to_3 = []
# df.apply(lambda x: nmdxtract(x['ID'], x['PTCposition'], gtfA, x['Variant_Classification']), axis=1)
#
# df['NMDA'] = NMD
# df['to_5A'] = to_5
# df['to_3A'] = to_3


df['NMDB'] = df.apply(lambda x: nmdxtract(x['ID'], x['PTCposition'], gtf, x['Variant_Classification']), axis=1)

# df['to_5B'] = to_5
# df['to_3B'] = to_3

# Save NMD-annotated df (NMD_Braun.csv)

df.to_csv('new_NMDB_NoUCSC_Braun2.0_HD_remix.csv',sep='\t',index = None)
