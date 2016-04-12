#!/usr/bin/env python

import xml.etree.ElementTree as ET

tree = ET.parse('../html/index.qhp')
root = tree.getroot()

for child in root:
  if child.tag == "filterSection":
    for grandchild in child:
      if grandchild.tag == "files":
        
        st = ET.Element("file")
        st.text = "SkirtLogoSmall-home.png"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/MathJax.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/extensions/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/extensions/TeX/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/extensions/HTML-CSS/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/fonts/HTML-CSS/TeX/eot/*.eot"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/fonts/HTML-CSS/TeX/otf/*.otf"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/element/mml/jax.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/element/mml/optable/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/input/TeX/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/autoload/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/AMS/Regular/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Caligraphic/Bold/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Fraktur/Bold/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Fraktur/Regular/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Greek/BoldItalic/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Main/Bold/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Main/Italic/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Main/Regular/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Math/BoldItalic/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/SansSerif/Bold/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/SansSerif/Italic/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/SansSerif/Regular/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Script/Regular/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/TeX/Typewriter/Regular/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/STIX/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/STIX/General/Bold/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/STIX/General/BoldItalic/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/STIX/General/Italic/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/STIX/General/Regular/*.js"
        grandchild.append(st)
        
        st = ET.Element("file")
        st.text = "mathjax/jax/output/HTML-CSS/fonts/STIX/IntegralsD/Regular/*.js"
        grandchild.append(st)
                
tree.write('../html/index_mathjax.qhp')
