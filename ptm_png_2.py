import cgi,cgitb
import os

import sys
import requests
import re
import json
import random
import datetime
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

cgitb.enable()

def get_protein(_l):
	url = 'http://rest.thegpm.org/1/protein/sequence/acc=%s' % (_l)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values = json.loads(r.text)
	except:
		return None
	return values[0]

def get_description(_l):
	a = _l
	if re.search(r'A.+\.\d+$',_l):
		a = re.sub(r'\.\d+$','',_l)
	url = 'http://rest.thegpm.org/1/protein/description/acc=%s' % (a)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values = json.loads(r.text)
	except:
		return None
	return values[0]

def make_ptm_csv(_l,_plength,_title,_protein,_file,_y,_s,_ltype):
	session = requests.session()
	seq = list(_protein)
	url = 'http://rest.thegpm.org/1/peptide/pf/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	values = {'acetyl':None,'phosphoryl':None,'ubiquitinyl':None}
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB') 
		return None
	try:
		values['phosphoryl'] = json.loads(r.text)
	except:
		print('JSON error: problem with phosphorylation data') 
		return None

	url = 'http://rest.thegpm.org/1/peptide/af/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB') 
		return None
	try:
		values['acetyl'] = json.loads(r.text)
	except:
		print('JSON error: problem with acetylation data') 
		return None

	url = 'http://rest.thegpm.org/1/peptide/uf/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB') 
		return None
	try:
		values['ubiquitinyl'] = json.loads(r.text)
	except:
		print('JSON error: problem with ubiquitinylation data') 
		return None
	#formulate a URL to request information about sumoylation for the protein identified by _l
	url = 'http://rest.thegpm.org/1/peptide/su/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values['K-sumoyl'] = json.loads(r.text)
	except:
		return None

	#formulate a URL to request information about R-dimethylation for the protein identified by _l
	url = 'http://rest.thegpm.org/1/peptide/di/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB') 
		return None
	try:
		values['dimethyl'] = json.loads(r.text)
	except:
		print('JSON error: problem with dimethylation data') 
		return None
	#formulate a URL to request information about O-linked glycosylation for the protein identified by _l
	url = 'http://rest.thegpm.org/1/peptide/ol/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB') 
		return None
	try:
		values['ST-glyco'] = json.loads(r.text)
	except:
		print('JSON error: problem with dimethylation data') 
		return None
	#formulate a URL to request information about R-citrullination for the protein identified by _l
	url = 'http://rest.thegpm.org/1/peptide/ct/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values['citrulline'] = json.loads(r.text)
	except:
		return None
	#formulate a URL to request information about KP-oxidation for the protein identified by _l
	url = 'http://rest.thegpm.org/1/peptide/ox/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB') 
		return None
	try:
		values['oxidation'] = json.loads(r.text)
	except:
		print('JSON error: problem with oxidation data') 
		return None
	a = 1;
	xs = {'citrulline':[],'acetyl':[],'S-phosphoryl':[],'T-phosphoryl':[],'Y-phosphoryl':[],'ubiquitinyl':[],'K-sumoyl':[],'R-dimethyl':[],'P-oxidation':[],'K-oxidation':[],'ST-glyco':[]}
	ys = {'citrulline':[],'acetyl':[],'S-phosphoryl':[],'T-phosphoryl':[],'Y-phosphoryl':[],'ubiquitinyl':[],'K-sumoyl':[],'R-dimethyl':[],'P-oxidation':[],'K-oxidation':[],'ST-glyco':[]}
	ts = {'acetyl':0,'phosphoryl':0,'oxidation':0,'ubiquitinyl':0,'K-sumoyl':0,'dimethyl':0,'citrulline':0,'ST-glyco':0}
	min_obs = 5
	lines = []
	ms = re.finditer(r'(?=N[^P][ST])',_protein)
	nl = ''
	if len([m for m in ms]) > 0:
		nl = '(<a class="bluesq" href="/_/nl_png/l=%s" target="_nl" title="Check for N-linked glycosylation">&#9632;</a>)' % (_l)
	phospho = {'S':0,'T':0,'Y':0}
	hydroxy = {'K':0,'P':0}
	oglyco = {'S':0,'T':0}
	kms = [_protein.find('M',1)]
	if kms[0] != -1:
		kms[0] = kms[0] + 1
		kms.append(kms[0]+1)
	for a in range(1,_plength+1):
		b = str(a)
		if b in values['acetyl']:
			if values['acetyl'][b] is None:
#				print('<script>document.getElementById("diagram").style="display: none;"</script>')
				print('<br />')
				print('<div><p>No PTMs detected for "%s".</p></div>' % (_l))
				return
		pre = ''
		if a - 4 >= 0:
			pre = seq[a-4].lower()
		elif a - 4 == -1:
			pre = '['
		if a - 3 >= 0:
			pre += seq[a-3].lower()
		elif a - 3 == -1:
			pre = '['
		if a - 2 >= 0:
			pre += seq[a-2].lower()
		elif a - 2 == -1:
			pre = '['
		post = ''
		if a < len(seq):
			post = seq[a].lower()
			if a+1 < len(seq):
				post += seq[a+1].lower()
			else:
				post += ']'
			if a+2 < len(seq):
				post += seq[a+2].lower()
			elif post.find(']') == -1:
				post += ']'
		elif a == len(seq):
			post = ']'
		if len(lines) % 2 != 0:
			line = '<tr><td>%i</td><td>%s</td>' % (a,'<i style="font-size: 10pt">%s&middot;</i>%s<i style="font-size: 10pt">&middot;%s</i>' % (pre,seq[a-1],post))
		else:
			line = '<tr class="alt"><td>%i</td><td>%s</td>' % (a,'<i style="font-size: 10pt">%s&middot;</i>%s<i style="font-size: 10pt">&middot;%s</i>' % (pre,seq[a-1],post))
		ok = False
		
		if(b in values['acetyl']):
			if values['acetyl'][b] >= min_obs and (seq[a-1] == 'K' or a < 4):
				xs['acetyl'].append(a)
				ys['acetyl'].append(values['acetyl'][b])
				ts['acetyl'] += 1
				line += '<td>%i</td>' % (values['acetyl'][b])
				ok = True
			elif(b in values['acetyl']) and a > 3 and a in kms:
				if values['acetyl'][b] >= min_obs:
					xs['acetyl'].append(a)
					ys['acetyl'].append(values['acetyl'][b])
					ts['acetyl'] += 1
					line += '<td>%i</td>' % (values['acetyl'][b])
					ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'
				
		if(b in values['ubiquitinyl']) and seq[a-1] == 'K':
			if values['ubiquitinyl'][b] >= min_obs:
				xs['ubiquitinyl'].append(a)
				ys['ubiquitinyl'].append(values['ubiquitinyl'][b])
				ts['ubiquitinyl'] += 1
				line += '<td>%i</td>' % (values['ubiquitinyl'][b])
				ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		if(b in values['K-sumoyl']):
			if values['K-sumoyl'][b] >= 2 and seq[a-1] == 'K':
				xs['K-sumoyl'].append(a)
				ys['K-sumoyl'].append(values['K-sumoyl'][b])
				ts['K-sumoyl'] += 1
				line += '<td>%i</td>' % (values['K-sumoyl'][b])
				ok = True
			else:
				line += '<td></td>'
				
		if(b in values['phosphoryl']):
			if values['phosphoryl'][b] >= min_obs:
				if seq[a-1] == 'S':
					xs['S-phosphoryl'].append(a)
					ys['S-phosphoryl'].append(values['phosphoryl'][b])
					ts['phosphoryl'] += 1
					phospho['S'] += 1
					line += '<td>%i</td>' % (values['phosphoryl'][b])
					ok = True
				elif seq[a-1] == 'T':
					xs['T-phosphoryl'].append(a)
					ys['T-phosphoryl'].append(values['phosphoryl'][b])
					ts['phosphoryl'] += 1
					phospho['T'] += 1
					line += '<td>%i</td>' % (values['phosphoryl'][b])
					ok = True
				elif seq[a-1] == 'Y':
					xs['Y-phosphoryl'].append(a)
					ys['Y-phosphoryl'].append(values['phosphoryl'][b])
					ts['phosphoryl'] += 1
					phospho['Y'] += 1
					line += '<td>%i</td>' % (values['phosphoryl'][b])
					ok = True
				else:
					line += '<td></td>'
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		if(b in values['dimethyl']):
			if values['dimethyl'][b] >= min_obs and seq[a-1] == 'R':
					xs['R-dimethyl'].append(a)
					ys['R-dimethyl'].append(values['dimethyl'][b])
					ts['dimethyl'] += 1
					line += '<td>%i</td>' % (values['dimethyl'][b])
					ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		if(b in values['oxidation']):
			if values['oxidation'][b] >= min_obs:
				if seq[a-1] == 'K':
					xs['K-oxidation'].append(a)
					ys['K-oxidation'].append(values['oxidation'][b])
					line += '<td>%i</td>' % (values['oxidation'][b])
					ts['oxidation'] += 1
					hydroxy['K'] += 1
					ok = True
				elif seq[a-1] == 'P':
					xs['P-oxidation'].append(a)
					ys['P-oxidation'].append(values['oxidation'][b])
					ts['oxidation'] += 1
					hydroxy['P'] += 1
					line += '<td>%i</td>' % (values['oxidation'][b])
					ok = True
				else:
					line += '<td></td>'
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'
			
		if(b in values['citrulline']):
			if values['citrulline'][b] >= min_obs and seq[a-1] == 'R':
					xs['citrulline'].append(a)
					ys['citrulline'].append(values['citrulline'][b])
					ts['citrulline'] += 1
					line += '<td>%i</td>' % (values['citrulline'][b])
					ok = True
			else:
				line += '<td></td>'
				
		else:
			line += '<td></td>'
		oglyco
		if(b in values['ST-glyco']):
			if values['ST-glyco'][b] >= min_obs and seq[a-1] == 'S':
				xs['ST-glyco'].append(a)
				ys['ST-glyco'].append(values['ST-glyco'][b])
				ts['ST-glyco'] += 1
				oglyco['S'] += 1
				line += '<td>%i</td>' % (values['ST-glyco'][b])
				ok = True
			elif values['ST-glyco'][b] >= min_obs and seq[a-1] == 'T':
				xs['ST-glyco'].append(a)
				ys['ST-glyco'].append(values['ST-glyco'][b])
				ts['ST-glyco'] += 1
				oglyco['T'] += 1
				line += '<td>%i</td>' % (values['ST-glyco'][b])
				ok = True
			else:
				line += '<td></td>'
				
		if ok:
			lines.append(line + '</tr>')

	if(_s == 'xkcd'):
		plt.xkcd()
	else:
		if len(_s) > 0:
			try:
				plt.style.use(_s)
				mpl.style.use(_s)
			except:
				mpl.style.use('seaborn-notebook')
		else:
			mpl.style.use('seaborn-notebook')
	plt.xlim(0,int(1.02*_plength))
	ms = 10
	ms8 = 8
	plt.plot(xs['acetyl'],ys['acetyl'],color=(0.25,0,1,.8),markersize=ms,marker='o',linestyle='None',label='acetyl')
	plt.plot(xs['S-phosphoryl'],ys['S-phosphoryl'],markersize=ms,color=(1,0,.25,.8),marker='v',linestyle='None',label='S-phos')
	plt.plot(xs['T-phosphoryl'],ys['T-phosphoryl'],markersize=ms,color=(1,0,.25,.8),marker='^',linestyle='None',label='T-phos')
	plt.plot(xs['Y-phosphoryl'],ys['Y-phosphoryl'],markersize=ms,color=(1,0,.25,.8),marker='o',linestyle='None',label='Y-phos')
	plt.plot(xs['ubiquitinyl'],ys['ubiquitinyl'],markersize=ms,color=(.1,.8,.1,.8),marker='v',linestyle='None',label='K-GG')
	plt.plot(xs['K-sumoyl'],ys['K-sumoyl'],markersize=ms,color=(.1,.8,.1,.8),marker='X',linestyle='None',label='K-sumo')
	plt.plot(xs['R-dimethyl'],ys['R-dimethyl'],markersize=ms,color=(.3,.3,.1,.8),marker='d',linestyle='None',label='R-dimet')
	plt.plot(xs['K-oxidation'],ys['K-oxidation'],markersize=ms,color=(.1,.8,.8,.8),marker='*',linestyle='None',label='K-oxy')
	plt.plot(xs['P-oxidation'],ys['P-oxidation'],markersize=ms,color=(.8,.1,.8,.8),marker='*',linestyle='None',label='P-oxy')
	plt.plot(xs['citrulline'],ys['citrulline'],markersize=ms,color=(.3,.3,.3,.8),marker='.',linestyle='None',label='R-citr')
	plt.plot(xs['ST-glyco'],ys['ST-glyco'],markersize=ms,markerfacecolor=(1,1,1,.8),markeredgewidth=1,markeredgecolor=(.6,.6,.6,.8),marker='s',linestyle='None',label='O-glyco')
	if _ltype == 'linear':
		plt.yscale('linear')
	else:
		plt.yscale('log')
	plt.ylabel('PSM (tabbs)')
	plt.xlabel('residue')
	plt.legend(loc='best')
	plt.grid(True, lw = 1, ls = '--', c = '.8')
#	plt.axvline(x=_plength,color=(.2,.2,.2,.5),linestyle='dotted',linewidth=1)
	plt.title(_title)
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	fig = plt.gcf()
	fig.set_size_inches(10, 5)
	cl = re.sub('[\|\:]','_',_file)
	plt.gca().get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	if _y is not None:
		print(_y)
		plt.ylim(1,_y)
	else:
		plt.ylim(1,None)
	ps = ax.get_ylim()
	if _ltype != 'linear':
		if ps[1] < 2000:
			ax.set_ylim([1,2000])
	desc = re.sub(r' \[',r'<br />[&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;',get_description(_l))
	desc = re.sub(r'[\[\]]',r'',desc)
	if len(lines) == 0:
		plt.ylim(0,1000)
	fig.savefig('../ptm_png/%s_ptms.png' % (cl), dpi=200, bbox_inches='tight')
	script = "<div id='diagram' class='pic'><center><img src='/ptm_png/%s_ptms.png' height='400' width='800' /></center></div>" % (cl)
	print(script)
	up = re.findall(r'UP\:(\w+)',desc)
	if len(up) and up[0] != 'NA':
		link = 'UP:<a href="https://uniprot.org/uniprot/%s" target="_blank">%s</a>; GlyGen: <a href="https://glygen.org/protein/%s-1#Glycosylation" target="_blank" title="GlyGen entry">%s</a>' % (up[0],up[0],up[0],up[0])
		desc = re.sub(r'UP\:\w+',link,desc)
	print('<p class="desc">%s %s  (<a href="http://gpmdb.thegpm.org/~/dblist_label/label=%s" title="main GPMDB page for this sequence" target="_blank">more</a>)</p>' % (re.sub(r'alt\:',r'<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;alt:',desc),nl,_l))
	create_table(lines,ts,_l,phospho,hydroxy,oglyco)
	return

def create_table(_lines,_ts,_l,_sty,_pk,_og):
	print('''
	<p class="con">This Modification-Abundance (M-A) diagram shows the number of times a residue has been observed with a particular PTM in a peptide-to-spectrum match in GPMDB, 
	as a function of the residue's position in the corresponding protein sequence. 
	Note that the Y-axis has a log scale. 
	The particular residues and observation numbers are detailed in the following table:</p>''')
	if len(_lines) == 0:
		print('<div><p style="margin: auto;border-collapse: collapse;">No PTMs detected</p></div>')
		return
	print('<div><table style="margin: auto;border-collapse: collapse;">')
	head = '''<tr class="heads"><td><span title="Location of modification in protein coordinates">pos</span></td>
	<td><span title="Amino acid residue with the modification">res</span></td>
	<td>&nbsp;&nbsp;&nbsp;<span title="PSMs with N-terminal or K-acetylation at a site (tabbs)">acetyl</span>&nbsp;&nbsp;&nbsp;</td>
	<td><span title="PSMs with  K-GG conjugation at a site (tabbs)">GG</span></td>
	<td><span title="PSMs with K-sumoylation at a site (tabbs)">sumo</span></td>
	<td><span title="PSMs with S/T/Y-phosphorylation at a site (tabbs)">phospho</span></td>
	<td><span title="PSMs with R-dimethylation at a site (tabbs)">R-dimeth</span></td>
	<td><span title="PSMs with P/K-oxidation at a site (tabbs)">P/K-oxy</span></td>
	<td><span title="PSMs with citrulline at a site (tabbs)">R-citr</span></td>
	<td><span title="PSMs with O-glycosylation at a site (tabbs)">O-glyco</span></td>
	</tr>'''
	print(head)
	for l in _lines:
		print(l)
	if len(_lines) > 10:
		print(head)
	pline = 'title="S=%i,T=%i,Y=%i"' % (_sty['S'],_sty['T'],_sty['Y'])
	kline = 'title="K=%i,P=%i"' % (_pk['K'],_pk['P'])
	gline = 'title="S+%i,T=%i"' % (_og['S'],_og['T'])
	print('''<tr class="tots"><td style="text-align: right">sites:</td>
	<td>%i</td>
	<td>%i</td>
	<td>%i</td>
	<td>%i</td>
	<td %s>%i</td>
	<td>%i</td>
	<td %s>%i</td>
	<td>%i</td>
	<td %s>%i</td>
	</tr>''' % (len(_lines),_ts['acetyl'],_ts['ubiquitinyl'],_ts['K-sumoyl'],pline,_ts['phosphoryl'],_ts['dimethyl'],kline,_ts['oxidation'],_ts['citrulline'],gline,_ts['ST-glyco']))
	print('</table></div>')
	print('''<div id="content"><ol>
	<li><b>pos:</b> location of the modification in protein coordinates;</li>
	<li><b>res:</b> amino acid residue with the modification and flanking residues;</li>
	<li><b>aceytl:</b> PSMs with N-terminal or K-acetylation at <i>pos</i> (tabbs);</li>
	<li><b>phospho:</b> PSMs with S/T/Y-phosphorylation at <i>pos</i> (tabbs);</li>
	<li><b>GG:</b> PSMs with K-GG conjugation at <i>pos</i> (tabbs);</li>
	<li><b>sumo:</b> PSMs with K-sumoylation at <i>pos</i> (tabbs);</li>
	<li><b>R-dimeth:</b> PSMs with R-dimethylation at <i>pos</i> (tabbs); &amp;</li>
	<li><b>P/K-oxy:</b> PSMs with 4-hydroxyproline or 5-hydroxylysine at <i>pos</i> (tabbs).</li>
	<li><b>R-citr:</b> PSMs with citrulline at <i>pos</i> (tabbs).</li>
	<li><b>O-glyco:</b> PSMs with O-linked glycosylation at <i>pos</i> (tabbs).</li>
	</ol></div>''')
	
form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
style = ''
try:
	style = form['s'].value
except:
	style = ''
	
ltype = 'log'
try:
	ltype = form['t'].value
except:
	ltype = 'log'

label = ''
try:
	label = form['l'].value
except:
	print('There must be a protein accession value specified')
	exit()
filename = label
title = 'ω-mod %s PTMs (%s)' % (label,ltype)
y_axis = None
protein = get_protein(label)
make_ptm_csv(label,len(protein),title,protein,filename,y_axis,style,ltype)
