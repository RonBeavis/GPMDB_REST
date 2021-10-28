#!c:/python3.7.4/python.exe

import cgi,cgitb
import re
import requests
import json
import os
import random
os.environ[ 'HOME' ] = 'c:/temp'
import sys
import subprocess
import time
import datetime
import re
import hashlib

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

def get_interpro(_l):
	ds = list()
	if not _l:
		return ds
	url = 'http://localhost/1/protein/domains/acc=%s' % (_l)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		ds = json.loads(r.text)
	except:
		return None
	return ds


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

def get_domains(_s,_hex):
	name = str(random.random())
	fasta = 'c:/temp/%s.fasta' % (name)
	f = open(fasta,'w')
	f.write('>%s no description\n' % (name))
	f.write(_s)
	f.close()
	param = ['c:/python3.7.4/Scripts/tmhmm.exe','-f',fasta,'-m','c:/python3.7.4/Scripts/TMHMM2.0.model']
	x = subprocess.Popen(param,
						stdin =subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True,
                        bufsize=0,cwd='c:/temp')
	while x.poll() == None:
		time.sleep(0.05)
	summary = 'c:/temp/%s.summary' % (name)
	vs = [re.sub(r' trans.+','',l.strip()).split(' ') for l in open(summary,'r') if l.find('transmembrane helix') != -1]
	os.remove(fasta)
	os.remove(summary)
	os.remove('c:/temp/%s.annotation' % (name))
	os.remove('c:/temp/%s.plot' % (name))
	ds = set()
	t = '';
	for i,v in enumerate(vs):
		if i == 0:
			t += '<div class="ex1">'
		elif i % 6 == 0:
			t += '<br />\n'
		for j in range(int(v[0])+1,int(v[1])+2):
			ds.add(j)
		t += '(%i-%i), ' % (int(v[0])+1,int(v[1])+1)
	if t:
		t = re.sub(r', $',r'',t)
		t += '</div>'
	return ds,t

def get_highlights(_h,_s):
	rs = set()
	if not _h:
		return rs
	h = _h.strip()
	vs = h.split(',')
	for v in vs:
		v = v.strip()
		if v.find('re:') == 0:
			p = re.compile(v[3:])
			for m in p.finditer(_s):
				for j in range(m.start()+1,m.end()+1):
					rs.add(j)
		elif v.find('-') != -1:
			gs = v.split('-')
			for j in range(int(gs[0]),int(gs[1])+1):
				rs.add(j)
		else:
			try:
				rs.add(int(v))
			except:
				try:
					p = re.compile('[%s]' % v)
					for m in p.finditer(_s):
						for j in range(m.start()+1,m.end()+1):
							rs.add(j)
				except:
					continue
	return rs

def get_marked(_h,_s):
	rs = set()
	if not _h:
		return rs
	h = _h.strip()
	vs = h.split(',')
	for v in vs:
		v = v.strip()
		if v.find('re:') == 0:
			p = re.compile(v[3:])
			for m in p.finditer(_s):
				for j in range(m.start()+1,m.end()+1):
					rs.add(j)
		elif v.find('-') != -1:
			gs = v.split('-')
			for j in range(int(gs[0]),int(gs[1])+1):
				rs.add(j)
		else:
			p = re.compile('[%s]' % v)
			for m in p.finditer(_s):
				for j in range(m.start()+1,m.end()+1):
					rs.add(j)
	return rs

def print_seq(_s,_m,_b,_l,_h,_fe):
	if len(_s) == 0:
		return
	m = hashlib.sha3_256()
	m.update(_s.encode())
	hexdigest = m.hexdigest()
	ds,text = get_domains(_s,hexdigest)
	hs = get_highlights(_h,_s)
	red = get_marked(_m,_s)
	blue = get_marked(_b,_s)
	interpro = get_interpro(_l)
	display = '<div class="ex1"><hr width="650" style="margin-left: -20px;"/>'
	mem = ''
	highlight = ''
	for i,r in enumerate(_s):
		if i != 0 and i % 50 == 0:
			display += '&nbsp;&nbsp;<span class="num">%i</span></div>\n<div class="ex1">' % (i)
		elif i % 10 == 0:
			display += ''
		mem = ''
		if i+1 in ds and _fe['grey']:
			mem = ' mem'
		highlight = ''
		if i+1 in hs and _fe['green']:
			highlight = ' highlight'			
		if i+1 in red and _fe['red']:
			display += '<span class="red%s%s" title="%s %i">%s</span>' % (highlight,mem,r,i+1,r.lower())
		elif i+1 in blue and _fe['blue']:
			display += '<span class="blue%s%s" title="%s %i">%s</span>' % (highlight,mem,r,i+1,r.lower())
		else:
			display += '<span class="unmarked%s%s" title="%s %i">%s</span>' % (highlight,mem,r,i+1,r.lower())
	display = re.sub(r'\<div class=\"ex1\"\>$',r'',display)
	display += '</div>\n'
	print(display)
	print('<div class="ex1"><hr width="650" style="margin-left: -20px;"/></div>\n')
	if text:
		print('<div class="ex1"><u>TM domains:</u></div>\n')
		print(text)
	if interpro:
		print('<div class="ex1"><u>Interpro domains:</u></div>\n')
		for v in interpro:
			print('<div class="ex1">(%i-%i) %s </div>' % (v['b'],v['e'],v['ldesc']))
	print('<div class="ex1"><u>SHA3 256:</u><br>\n%s\n</div>\n' % (hexdigest))
	return

def print_form(_s,_m,_b,_l,_h,_fe):
	but = '<input type="submit" class="button" value="&#8635;" title="refresh display" />'
	print('<div class="ex1">')
	print('<form style="display: inline;" name="seq_form" action="/thegpm-cgi/seq.py" METHOD="POST" ENCTYPE="multipart/form-data">')
	print('<hr width="650" style="margin-left: -20px;"/>')
	print('<input type="hidden" value="no" name="red" />')
	print('<input type="hidden" value="no" name="blue" />')
	print('<input type="hidden" value="no" name="green" />')
	print('<input type="hidden" value="no" name="grey" />')
	if _fe['grey']:
		print('<span class="mem"><input type="checkbox" id="grey_box" name="grey" value="yes" CHECKED/>&nbsp;TM domains</span>&nbsp;%s<br />' % (but))
	else:
		print('<span class="mem"><input type="checkbox" id="grey_box" name="grey" value="yes" />&nbsp;TM domains </span>&nbsp;%s<br />' % (but))
	if _fe['red']:
		print('<span class="red"><input type="checkbox" id="red_box" name="red" value="yes" CHECKED/>&nbsp;&nbsp;residues:</span>&nbsp;<input id="red" name="m" size="20" value="%s" placeholder="KR" />&nbsp;%s<br/>' % (_m,but))
	else:
		print('<span class="red"><input type="checkbox" id="red_box" name="red" value="yes"/>&nbsp;&nbsp;residues:</span>&nbsp;<input id="red" name="m" size="20" value="%s" placeholder="KR" />&nbsp;%s<br/>' % (_m,but))
	if _fe['blue']:
		print('<span class="blue"><input type="checkbox" id="blue_box" name="blue" value="yes" CHECKED/>&nbsp;&nbsp;residues:</span>&nbsp;<input id="blue" name="b" size="20" value="%s" placeholder="ED" />&nbsp;%s<br/>' % (_b,but))
	else:
		print('<span class="blue"><input type="checkbox" id="blue_box" name="blue" value="yes" />&nbsp;&nbsp;residues:</span>:&nbsp;<input id="blue" name="b" size="20" value="%s" placeholder="ED" />&nbsp;%s<br/>' % (_b,but))
	if _fe['green']:
		print('<span class="highlight"><input type="checkbox" id="green_box" name="green" value="yes" CHECKED/>&nbsp;&nbsp;&nbsp;&nbsp;ranges:</span>&nbsp;<input id="blue" name="h" size="20" value="%s" placeholder="1-20,25 or re:N[^P][ST]" />&nbsp;%s<br/>' % (_h,but))
	else:
		print('<span class="highlight"><input type="checkbox" id="green_box" name="green" value="yes" />&nbsp;&nbsp;&nbsp;&nbsp;ranges:</span>&nbsp;<input id="blue" name="h" size="20" value="%s" placeholder="1-20,25 or re:N[^P][ST]" />&nbsp;%s<br/>' % (_h,but))
	print('<span class="accession">&nbsp;&nbsp;&nbsp;accession:</span>&nbsp;<input id="label" name="l" size="20" value="%s" onChange="clearSeq();" placeholder="ENSP00000242786" />&nbsp;%s<br/>' % (_l,but))
	print('<span class="accession">&nbsp;&nbsp;&nbsp;&nbsp;sequence:</span><br/><textarea rows="10" cols="50" id="seq" name="s" placeholder="protein sequence">%s</textarea>&nbsp;%s' % (_s,but))
	print('</form></div>\n')
	
cgitb.enable()
form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
seq = ''
form_entries = dict()
try:
	seq = form['s'].value.upper()
except:
	seq = ''
seq = re.sub(r'[^A-Z]+',r'',seq)
mark = ''
try:
	mark = form['m'].value
except:
	mark = ''
form_entries['red'] = True
blue = ''
try:
	blue = form['b'].value
except:
	blue = ''
form_entries['blue'] = True
label = ''
try:
	label = form['l'].value.upper()
except:
	label = ''
highlight = ''
try:
	highlight = form['h'].value
except:
	highlight = ''
boxes = ['grey','red','green','blue']
for b in boxes:
	try:
		if form.getvalue(b) == 'yes':
			form_entries[b] = True
		else:
			form_entries[b] = False
	except:
		form_entries[b] = True
if seq:
	print_seq(seq,mark,blue,label,highlight,form_entries)
if not seq and label:
	seq = get_protein(label)
	print_seq(seq,mark,blue,label,highlight,form_entries)
print_form(seq,mark,blue,label,highlight,form_entries)

