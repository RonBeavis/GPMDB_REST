#!c:/python3.7.4/python.exe

import cgi,cgitb
import sys
import requests
import re
import json
import datetime

def get_peptides(_e):
	url = 'http://gpmdb.thegpm.org/protein/model/%s&excel=1' % (_e)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None

	text = re.sub('\r\n','\n',r.text)
	return text.splitlines()

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
	url = 'http://rest.thegpm.org/1/protein/description/acc=%s' % (_l)
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

def print_top(_l,_l2):
	desc = "Sequence overlap %s ∩ %s" % (_l,_l2)

	print('''<!DOCTYPE html>
		<html lang="en" class="no-js">
		<head>
		<meta http-equiv="X-UA-Compatible" content="IE=edge">
		<meta charset="utf-8">
		<title>Sequence overlap/intersection display</title>
		<meta name="viewport" content="width=device-width,initial-scale=1" />
		<meta name="robots" content="index,nofollow,noarchive">''')
	print('''
		<meta property="og:locale" content="en_EN" />
		<meta property="og:type" content="website" />
		<meta property="og:title" content="GPMDB Sequence overlap/intersection display" />
		<meta property="og:description" content="%s" />
		<meta property="og:url" content="https://gpmdb.thegpm.org" />
		<meta property="og:image:width" content="800" />
		<meta property="og:image:height" content="400" />
		<meta property="og:image" content="https://gpmdb.thegpm.org/pics/gpmdb_sq.png" />
		<meta property="og:image:secure_url" content="https://gpmdb.thegpm.org/pics/gpmdb_sq.png" />
		''' % (desc))
	v = re.sub(r'[\|\:]',r'_',_l)
	print('''
		<meta name="twitter:url" content="https://gpmdb.thegpm.org/_/overlap/l=%s&l2=%s">
		<meta name="twitter:domain" content="gpmdb.thegpm.org">
		<meta name="twitter:card" content="summary_large_image" />
		<meta name="twitter:site" content="@norsivaeb" />
		<meta name="twitter:description" content="%s" />
		<meta name="twitter:title" content="GPMDB Sequence overlap display - %s" />
		<meta name="twitter:image" content="https://gpmdb.thegpm.org/pics/gpmdb_sq.png" />
		'''  % (_l,_l2,desc,desc))
	print('''
		<style media="screen" type="text/css">
		@font-face	{
			font-family: 'Anonymous Pro';
			font-style: normal;
			font-weight: 400;
			src: local('Anonymous Pro'), local('Anonymous Pro-Regular'), url('https://gpmdb.thegpm.org/fonts/AnonymousPro-Regular.ttf');
			format('ttf');
		}
		@font-face	{
			font-family: 'Anonymous Pro';
			font-style: normal;
			font-weight: 700;
			src: local('Anonymous Pro-Bold'), url('https://gpmdb.thegpm.org/fonts/AnonymousPro-Bold.ttf');
			format('ttf');
		}
		@font-face	{
			font-family: 'Anonymous Pro';
			font-style: italic;
			font-weight: 400;
			src: local('Anonymous Pro-Italic'), url('https://gpmdb.thegpm.org/fonts/AnonymousPro-Italic.ttf');
			format('ttf');
		}
		body {
			color: #000000;
			background-color: #FFFFFF;
			font-weight: normal;
			font-family: "Anonymous Pro",serif;
			font-size: 13pt;
			margin: auto;
		}
		.cdiv	{
			  display: table;
			  margin: 0 auto;
		}
		.num	{
			color: grey;
			font-size: 11pt;
		}
		.red	{
			background-color: #ff6666;
			color: white;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.blue	{
			background-color: #6666ff;
			color: white;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.accession	{
			background-color: #996633;
			color: white;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.unmarked	{
			color: grey;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.mem	{
			background-color: #aaaaaa;
			color: #FFFFFF;
		}
		.highlight	{
			background-color: #00cc99;
			color: white;
		}
		div.ex1	{
			margin: 3px 3px 3px 3px;
		}
		.button {
		  font-weight: normal;
		  position: relative;
		  background-color: #4CAF50;
		  border: none;
		  font-size: 16px;
		  color: #FFFFFF;
		  padding: 3px;
		  width: 30px;
		  text-align: center;
		  -webkit-transition-duration: 0.4s; /* Safari */
		  transition-duration: 0.4s;
		  text-decoration: none;
		  overflow: hidden;
		  cursor: pointer;
		  margin-bottom: 8px;
		  margin-top: 8px;

		}

		.button:after {
			font-weight: normal;
		  content: "";
		  background: #f1f1f1;
		  display: block;
		  position: absolute;
		  padding-top: 300%;
		  padding-left: 10%;
		  margin-left: -20px !important;
		  margin-top: -120%;
		  opacity: 0;
		  transition: all 0.8s;
		}

		.button:active:after {
			font-weight: normal;
		  padding: 0;
		  margin: 0;
		  opacity: 1;
		  transition: 0s
		}
	</style>
	</head>''')
	print('''\n<body>
	<div id="main"><div id="main_body" class="cdiv">\n''')
	return

def print_bottom():
	print('<p id="copyright">%s GPMDB&nbsp;<img style="vertical-align: middle;" src="/pics/the_gpm_db_140x105.png" height="20" /></p>' % (datetime.datetime.now()))

	t = '''</div></div></body></html>\n'''
	print(t)
	return

def print_form(_l1,_l2):
	but = '<input type="submit" class="button" value="&#8635;" title="refresh display" />'
	print('<div class="ex1">')
	print('<form style="display: inline;" name="seq_form" action="/_/∩/" METHOD="GET">')
	print('<input id="red" name="l" size="20" value="%s" placeholder="ENSP0...." /> overlaps with <input name="l2" size="20" value="%s" placeholder="ENSP0...." />' % (_l1,_l2))
	print(but)
	print('</form>')

cgitb.enable()

form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
try:
	label1 = form['l'].value
	label2 = form['l2'].value
	print_top(label1,label2)
except:
	print_top('','')
	print_form('','')
	print('</body></html>\n')
	exit()

ls = get_peptides(label1)
protein_1 = {}
for l in ls:
	v = l.strip()
	if v.find('Sequence') != -1:
		continue
	vs = v.split('\t')
	if len(vs) < 5:
		continue
	if vs[5] not in protein_1:
		protein_1[vs[5]] = [(int(vs[0]),int(vs[1]))]
	else:
		protein_1[vs[5]].append((int(vs[0]),int(vs[1])))

protein_2 = {}
ls = get_peptides(label2)

mset = set()
for l in ls:
	v = l.strip()
	if v.find('Sequence') != -1:
		continue
	vs = v.split('\t')
	if len(vs) < 5:
		continue
	if vs[5] not in protein_1:
		continue
	if vs[5] not in protein_2:
		protein_2[vs[5]] = [(int(vs[0]),int(vs[1]))]
	else:
		protein_2[vs[5]].append((int(vs[0]),int(vs[1])))

for p in protein_2:
	ls = protein_2[p]
	for l in ls:
		for i in range(l[0],l[1]+1):
			mset.add(i)

seq = get_protein(label2)
des = re.sub('\[.+',r'',get_description(label2))
print('<p><a href="/protein/model/%s" target="_blank">%s</a>&mdash;%s <br/>' % (label2,label2,des))
print('overlap with observed peptides in <br/>')
des = re.sub('\[.+',r'',get_description(label1))
print('<a href="/protein/model/%s" target="_blank">%s</a>&mdash;%s </p>' % (label1,label1,des))
print('<p><span class="red">overlapping residues</span>: %i<br/>total residues: %i<br/><span class="red">overlap</span>: %.1f%%</p>' % (len(mset),len(seq),100*(len(mset)/len(seq))))
#print('<pre>')
display = '<div class="ex1">\n'
for i,s in enumerate(seq):
	if i != 0 and i % 50 == 0:
			display += '&nbsp;&nbsp;<span class="num">%i</span></div>\n<div class="ex1">' % (i)
	if i+1 in mset:
		display += '<span class="red" title="%s %i">%s</span>' % (s,i+1,s)
	else:
		display += '<span class="unmarked" title="%s %i">%s</span>' % (s,i+1,s)
display += '</div>\n'
print(display)
print_form(label1,label2)
print_bottom()

