#!c:/python3.7.4/python.exe

import cgi,cgitb
import os
os.environ[ 'HOME' ] = 'c:/temp'

import sys
import requests
import re
import json
import random
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

cgitb.enable()

def error(_x,_y):
	dx = np.sqrt(_x)
	dy = np.sqrt(_y)
	return np.sqrt((dx/_x)**2 + (dy/_y)**2 - 2*(dx/_x)*(dy/_y))
#
# get peptide PSMs for protein _l
#

def get_peptides(_l):
	a = _l
	if re.search(r'A.+\.\d+$',_l):
		a = re.sub(r'\.\d+$','',_l)
	url = 'http://localhost/protein/model/%s&excel=1' % (a)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None

	text = re.sub('\r\n','\n',r.text)
	return text.splitlines()
#
# get protein _l sequence
#

def get_protein(_l):
	url = 'http://localhost/1/protein/sequence/acc=%s' % (_l)
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
#
# get protein _l description
#
def get_description(_l):
	url = 'http://localhost/1/protein/description/acc=%s' % (_l)
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
#
# test deamidated residue for sequence motifs
#
def test_pos(_res,_protein):
	if _protein[_res+1] == 'G':
		return 50
	if _protein[_res+1] == 'A':
		return 25
	if _protein[_res-1] == 'K' or _protein[_res-1] == 'R':
		return 20
	if _protein[_res-1] == 'N' or _protein[_res+1] == 'N':
		return 20
	if _protein[_res-1] == 'Q' or _protein[_res+1] == 'Q':
		return 20
	return 15

def make_ptm_csv(_l,_plength,_title,_protein,_xs,_ys,_min,_legend):
	use_ylim = False
	notes_min = 0
	session = requests.session()
	seq = list(_protein)
	values = {'deamidation':None}
	#formulate a URL to request information about NQ-deamidation for the protein identified by _l
	url = 'http://localhost/1/peptide/nq/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values['deamidation'] = json.loads(r.text)
	except:
		return None
	a = 1;
	#xs and ys contain the x,y arrays for the scatter plots
	xs = {'N-deamidation':[],'N':[]}
	ys = {'N-deamidation':[],'N':[]}
	min_obs = 2
	#create x,y arrays for plot
	for a in range(1,_plength+1):
		b = str(a)
		if(b in values['deamidation']):
			if values['deamidation'][b] >= min_obs and a < len(seq) - 2:
				if seq[a-1] == 'N' and (seq[a+1] == 'S' or seq[a+1] == 'T') and seq[a] != 'P':
					xs['N-deamidation'].append(a)
					ys['N-deamidation'].append(values['deamidation'][b]) 
			elif seq[a-1] == 'N' and _ys[a-1] > 0:
				xs['N'].append(a)
				ys['N'].append(0)
	notes = []
	tsv = {}
	isv = {}
	for i,y in enumerate(ys['N-deamidation']):
		yst = ys['N-deamidation'][i]
		ystm = _ys[xs['N-deamidation'][i]-1]
		try:
			ys['N-deamidation'][i] /= _ys[xs['N-deamidation'][i]-1]/100.0
		except:
			notes.append('exception')
			continue
		if ys['N-deamidation'][i] >= 100.0:
			ys['N-deamidation'][i] = 99.0
		s = xs['N-deamidation'][i]
		if s - 2 < 0:
			p = '[' + _protein[s-1] + _protein[s].lower()
		elif s > len(_protein) - 1:
			p = _protein[s-2].lower() + _protein[s-1] + ']'
		else:
			p = _protein[s-2].lower() + _protein[s-1] + _protein[s].lower() + _protein[s+1].lower()
		tsv[xs['N-deamidation'][i]] = '%i\t%s\t%.2f\t%i\t%i\t%.2f' % (xs['N-deamidation'][i],p,ys['N-deamidation'][i],yst,ystm,error(yst,ystm))
		isv[xs['N-deamidation'][i]] = ys['N-deamidation'][i]
		if s - 2 < 0:
			notes.append(('%sN%i%s'%('[',xs['N-deamidation'][i],_protein[s].lower()+_protein[s+1].lower()),xs['N-deamidation'][i],ys['N-deamidation'][i]))
		elif s > len(_protein) - 1:
			notes.append(('%sN%i%s'%(_protein[s-2].lower(),xs['N-deamidation'][i],']'),xs['N-deamidation'][i],ys['N-deamidation'][i]))
		else:
			notes.append(('%sN%i%s'%(_protein[s-2].lower(),xs['N-deamidation'][i],_protein[s].lower()+_protein[s+1].lower()),xs['N-deamidation'][i],ys['N-deamidation'][i]))
	mpl.style.use('seaborn-notebook')
	plt.xlim(0,int(1.02*_plength))
	max_occ = 0
	if len(ys['N-deamidation']) > 0:
		max_occ = max(ys['N-deamidation'])
		
	if use_ylim and max_occ < 9.5:
		plt.ylim(0,10.1)
	ms = 10
	tlist = ys['N-deamidation'] + ys['N']
	ave = 0.0
	if len(tlist) > 0:
		thelist = []
		for t in tlist:
			if t > 2.0:
				continue
			thelist.append(t)
		if len(thelist):
			ave = sum(thelist)/len(thelist)
	#load all x,y information into a plot
	xsf = []
	ysf = []
	ns = []
	for i,f in enumerate(xs['N-deamidation']):
		lim = test_pos(f-1,protein)
		if ys['N-deamidation'][i] >= lim and _ys[f-1] > 40:
			ysf.append(ys['N-deamidation'][i])
			xsf.append(f)
			ns.append(notes[i])
	notes = ns
	plt.plot(xsf,ysf,markersize=ms,color=(0.05,.05,.9,.8),marker='s',linestyle='None',label='N-linked')
	if len(xsf) == 0:
		plt.ylim(0,100.0)
	#set up the required graph
	plt.legend(loc=_legend)
	plt.yscale('linear')
	plt.ylabel('AI-ND score')
	plt.xlabel('residue')
	plt.grid(True, lw = 1, ls = '--', c = '.8')
	plt.title(_title)
	ax = plt.gca()
	xoff = 0.015 * len(protein)
	for n in notes:
		ax.annotate(n[0], (n[1]+xoff, n[2]))
	box = ax.get_position()
	ax.set_ylim([0,100])
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))


	fig = plt.gcf()
	fig.set_size_inches(10, 5)
	plt.gca().get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	#store graph in png file
	desc = re.sub(r' \[',r'<br />[&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;',get_description(_l))
	desc = re.sub(r'[\[\]]',r'',desc)
	up = re.findall(r'UP\:(\w+)',desc)
	if False and len(xsf) == 0:
		ms = re.finditer(r'(?=N[^P][ST])',protein)
		rs = [m for m in ms]
		if len(rs) == 1:
			print('<p class="con">No N-linked glycosylation sites observed (0 of %i <i>N[^P][ST]</i> site).</p>' % (len(rs)))
		else:
			print('<p class="con">No N-linked glycosylation sites observed (0 of %i <i>N[^P][ST]</i> sites).</p>' % (len(rs)))
		return
	if len(up) and up[0] != 'NA':
		link = 'UP:<a href="https://uniprot.org/uniprot/%s" target="_blank">%s</a>; GlyGen: <a href="https://glygen.org/protein/%s-1#Glycosylation" target="_blank" title="GlyGen entry">%s</a>' % (up[0],up[0],up[0],up[0])
		desc = re.sub(r'UP\:\w+',link,desc)
	fig.savefig('../ptm_png/%s_nl.png' % (re.sub(r'[\|\:]',r'_',_l)), dpi=200, bbox_inches='tight')
	script = "<img src='/ptm_png/%s_nl.png' height='400' width='800' />" % (re.sub(r'[\|\:]',r'_',_l))
	print("<div id='diagram' class='pic'><center>%s</center></div>" % (script))
	nl = '(<a class="bluesq" href="/_/ptm_png/l=%s" target="_ptm" title="Check for common PTMs">&#9679;</a>)' % (_l)
	desc
	print('<p class="desc">%s %s</p>' % (re.sub(r'alt\:',r'<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;alt: ',desc),nl))
	lines = []
	for i,x in enumerate(xsf):
		if len(lines) % 2 == 0:
			line = '<tr><td>%i</td><td><i style="font-size: 10pt">' % (x)
		else:
			line = '<tr class="alt"><td>%i</td><td><i style="font-size: 10pt">' % (x)
		if x-4 >= 0:
			line += '%s' % (protein[x-4].lower())
		else:
			line += '%s' % ('[')
		if x-3 >= 0:
			line += '%s' % (protein[x-3].lower())
		else:
			line += '%s' % ('[')
		if x-2 >= 0:
			line += '%s' % (protein[x-2].lower())
		else:
			line += '%s' % ('[')
		try:
			line += '&middot;</i>%s<i style="font-size: 10pt">&middot;%s%s%s</i></td><td>%.1f</td><td>%.0f</td></tr>' % (protein[x-1],protein[x].lower(),protein[x+1].lower(),protein[x+2].lower(),ysf[i],_ys[x-1])
		except:
			line += '&middot;</i>%s<i style="font-size: 10pt">&middot;%s%s%s</i></td><td>%.1f</td><td>%.0f</td></tr>' % (protein[x-1],protein[x].lower(),protein[x+1].lower(),']',ysf[i],_ys[x-1])

		lines.append(line)
	create_table(lines,protein)
	return

def create_table(_lines,_protein):
	ms = re.finditer(r'(?=N[^P][ST])',_protein)
	rs = [m for m in ms]
	if len(_lines) == 0:
		if len(rs) == 1:
			print('<p class="con">No N-linked glycosylation sites observed (0 of %i <i>N[^P][ST]</i> site).</p>' % (len(rs)))
		else:
			print('<p class="con">No N-linked glycosylation sites observed (0 of %i <i>N[^P][ST]</i> sites).</p>' % (len(rs)))
		return
	print('''
	<p class="con">This N-linked glysolyation site identification diagram shows the residues in a protein most likely to have this modification, 
	based on the data in GPMDB. The y-axis value indicates the confidence of the assignment.</p>''')
	print('<div><table style="margin: auto;border-collapse: collapse;">')
	print('''<tr class="tots"><td><span title="Location of modification in protein coordinates">pos</span></td>
	<td><span title="Amino acid residue with the modification (upper case) with flanking residues (lower case)">res</span></td>
	<td><span title="site assignment score: higher is better">AI-ND</span></td>
	<td><span title="PSMs used to identify this site (tabbs)">PSMs (Ta)</span></td>
	</tr>''')
	for l in _lines:
		print(l)
	print('''<tr class="tots"><td style="text-align: right">sites:</td>
	<td>%i (of %i)</td><td></td><td></td></tr>''' % (len(_lines),len(rs)))
	print('</table></div>')
	print('''<div id="content"><ol>
	<li><b>pos:</b> location of the modification in protein coordinates;</li>
	<li><b>res:</b> amino acid residue with the modification (upper case) with flanking residues (lower case);</li>
	<li><b>AI-ND:</b> assignment confidence (higher is better, ranging from 5 &ndash;100); and</li>
	<li><b>PSMs:</b> the number of peptide-to-spectrum assignments used to identify a site (tabbs).</li>
	</ol></div>''')
	
#
# retrieve the number of times each residue has been observed
#

def get_residues(_plength,_lines):
	start = {}
	end = {}
	res = {}
	max = 0
	for line in _lines[1:]:
		vs = line.split('\t')
		vs[0] = int(vs[0])
		vs[1] = int(vs[1])
		vs[2] = int(vs[2])
		if vs[0] in start:
			start[vs[0]] = start[vs[0]] + vs[2]
		else:
			start[vs[0]] = vs[2]
		if vs[1] in end:
			end[vs[1]] = end[vs[1]] + vs[2]
		else:
			end[vs[1]] = vs[2]
		if vs[1] > max:
			max = vs[1]
		for a in range(vs[0],vs[1]+1):
			if a in res:
				res[a] = res[a] + vs[2]
			else:
				res[a] = vs[2]

	max = _plength
	max_res = 0
	for a in res:
		if res[a] > max_res:
			max_res = res[a]
	a = 1
	xs = []
	ys = []
	for a in range(1,max+1):
		xs.append(a)
		if a in res:
			ys.append(res[a])
		else:
			ys.append(0)
	return (xs,ys)

form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
try:
	label = form['l'].value
except:
	print('There must be a protein accession value specified')
	exit()

filename = label
title = '%s N-linked sites' % (label)
y_axis = None

protein = get_protein(label)
min = 5
ls = get_peptides(label)
(xs,ys) = get_residues(len(protein),ls)
legend = 'upper left'
make_ptm_csv(label,len(protein),title,protein,xs,ys,min,legend)


