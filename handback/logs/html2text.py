"""Plain-text dump of an HTML file (drops <script>/<style> and HTML comments). Used only for quoting web pages."""
import sys, re, html
from html.parser import HTMLParser
class P(HTMLParser):
    def __init__(s):
        super().__init__(convert_charrefs=True); s.out=[]; s.skip=0
    def handle_starttag(s,t,a):
        if t in("script","style"): s.skip+=1
        if t in("br","p","div","tr","li","h1","h2","h3","h4","h5","h6","table","ul","ol","td","th"): s.out.append("\n")
        if t=="a":
            h=dict(a).get("href")
            if h: s.out.append(f" [link:{h}] ")
    def handle_endtag(s,t):
        if t in("script","style"): s.skip-=1
        if t in("p","div","tr","li","h1","h2","h3","h4","h5","h6","table"): s.out.append("\n")
        if t in("td","th"): s.out.append(" | ")
    def handle_data(s,d):
        if not s.skip: s.out.append(d)
raw=open(sys.argv[1],"rb").read().decode(sys.argv[2] if len(sys.argv)>2 else "utf-8",errors="replace")
p=P(); p.feed(raw)
t="".join(p.out)
t=re.sub(r"[ \t\r\xa0]+"," ",t); t=re.sub(r" *\n *","\n",t); t=re.sub(r"\n{3,}","\n\n",t)
print(t.strip())
