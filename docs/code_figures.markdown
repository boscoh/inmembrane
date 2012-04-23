***

<pre>
lcl|AE004092.1_cdsid_AAK34689.1   PSE-Lipoprotein     lipop;signalp                                 AE004092.1_cdsid_AAK34689.1 [gene=lmb]
lcl|AE004092.1_cdsid_AAK34690.1   PSE-Cellwall        hmm(Gram_pos_anchor,LPxTG);signalp;tmhmm(2)   AE004092.1_cdsid_AAK34690.1 [gene=SPy_2009]
lcl|AE004092.1_cdsid_AAK34691.1   PSE-Cellwall        hmm(LPxTG_PS50847);signalp                    AE004092.1_cdsid_AAK34691.1 [gene=scpA]
lcl|AE004092.1_cdsid_AAK34692.1   CYTOPLASM(non-PSE)  .                                             AE004092.1_cdsid_AAK34692.1 [gene=SPy_2013]
lcl|AE004092.1_cdsid_AAK34693.1   SECRETED            signalp;tmhmm(1)                              AE004092.1_cdsid_AAK34693.1 [gene=sic]
lcl|AE004092.1_cdsid_AAK34694.1   PSE-Cellwall        hmm(Gram_pos_anchor,LPxTG);signalp;tmhmm(1)   AE004092.1_cdsid_AAK34694.1 [gene=emm1]
lcl|AE004092.1_cdsid_AAK34695.1   CYTOPLASM(non-PSE)  .                                             AE004092.1_cdsid_AAK34695.1 [gene=mga]
lcl|AE004092.1_cdsid_AAK34696.1   MEMBRANE(non-PSE)   tmhmm(2)                                      AE004092.1_cdsid_AAK34696.1 [gene=SPy_2023]
lcl|AE004092.1_cdsid_AAK34697.1   SECRETED            signalp;tmhmm(1)                              AE004092.1_cdsid_AAK34697.1 [gene=isp]
lcl|AE004092.1_cdsid_AAK34698.1   PSE-Membrane        tmhmm(1)                                      AE004092.1_cdsid_AAK34698.1 [gene=SPy_2026]
lcl|AE004092.1_cdsid_AAK34699.1   CYTOPLASM(non-PSE)  .                                             AE004092.1_cdsid_AAK34699.1 [gene=SPy_2027]
</pre>

> Figure 2. Example of output

****


******

```python
    if is_hmm_profile_match:
      category =  "PSE-Cellwall"
    elif has_tm_helix(protein):
      if has_surface_exposed_loop(protein):
        category = "PSE-Membrane"
      else:
        category = "MEMBRANE(non-PSE)"
    else:
      if is_lipop:
        # whole protein considered outer terminal loop
        if sequence_length(protein) < terminal_exposed_loop_min:
          category = "LIPOPROTEIN(non-PSE)"
        else:
          category = "PSE-Lipoprotein"
      elif is_signalp:
        category = "SECRETED"
      else:
        category = "CYTOPLASM(non-PSE)"
```

> Figure 3. Main logic classifying subcellular localization and potential surface exposure for Gram-positive protein sequences, expressed in Python code. This algorithm was adapted from _SurfG+_. The function `has_surface_exposed_loop` evaluates whether the extracellular loops are sufficifiently long to be exposed out of the peptidoglycan layer. The rule adapted from _SurfG+_ is a minimum length of 50 amino acids for terminal loops, and 100 amino acids for internal loops.

*******


***********

```python
    def annotate(params, proteins):
      for seqid in proteins:
        proteins[seqid]['is_signalp'] = False
        proteins[seqid]['signalp_cleave_position'] = None

      signalp4_out = 'signalp.out'
      cmd = '%(signalp4_bin)s -t %(signalp4_organism)s  %(fasta)s' % \
                 params
      helpers.run(cmd, signalp4_out)
      
      for line in open(signalp4_out):
        if line.startswith("#"):
          continue
        words = line.split()
        seqid = helpers.parse_fasta_header(words[0])[0]
        if (words[9] == "Y"):
          proteins[seqid]['is_signalp'] = True
          proteins[seqid]['signalp_cleave_position'] = int(words[4])

      return proteins
```    

> Figure 4. Example of parsing code in the signalp4 plugin. The entire function responsible for processing _SignalP_ output. `helpers` is an _inmembrane_ module with utility functions. 

***********

***********

```python
>>> from twill.commands import *                                                 (1)
>>> go("http://services.cbu.uib.no/tools/bomp/")                                 (2)
==> at http://services.cbu.uib.no/tools/bomp/
'http://services.cbu.uib.no/tools/bomp/'
>>> showforms()                                                                  (3)

Form #1
## ## __Name__________________ __Type___ __ID________ __Value__________________
1     seqs                     textarea  (None)        
2     useblast                 checkbox  blast        [] of ['on'] 
3     evalue                   select    (None)       ['0.0000000001'] of ..
4     queryfile                file      (None)       None 
5  1  SUBMIT                   submit    (None)       Submit Search 
6     None                     reset     (None)       None 

..
>>> formfile("1", "queryfile", "/home/perry/my_sequences.fasta")                 (4)

>>> submit()                                                                     (5)

>>> links = showlinks()                                                          (6)
Links:
0. UiB ==> http://www.uib.no
..snip..
10. services@cbu.uib.no ==> mailto:services@cbu.uib.no
11. Click to check output status (new window) ==> viewOutput?id=99941690
12.  ==> http://www.unifob.uib.no/
13.  ==> http://www.uib.no/
..
>>> go("viewOutput?id=99941690")                                                 (7)
==> at http://services.cbu.uib.no/tools/bomp/viewOutput?id=99941690
'http://services.cbu.uib.no/tools/bomp/viewOutput?id=99941690'
>>> out = show()                                                                 (8)
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" 
 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
..
```

> Figure 5. An example of interfacing with the BOMP ß-barrel outer membrane protein predictor (Berven et al, 2004) web site using _twill_ on the Python interactive commandline. _twill_ essentially behaves like a headless web-browser. Lines with `>>>` denote inputs to the Python interactive command line, while other lines are output from _twill_ (1) First the appropriate commands from the _twill_ library are imported. (2) We navigate to the BOMP website, which silently downloads the HTML page and (3) show a summary of the forms on that page, including field names and input types. (4) We then use the `formfile` function to associate a local file with the `queryfile` FILE input field. Calling `submit()` (5) is equivalent to clicking the SUBMIT button defined in the form. After a short delay, an intermediate page is returned, and we can list the hyperlinks on this page using (6) showlinks(), and assign them to a variable (`links`, a Python list). We can then navigate to the appropriate result page (7) and assign the HTML text of this page to a variable (`out`) (8) for downstream parsing using BeautifulSoup. This type of interactive exploration can be easily expanded into an _inmembrane_ plugin to programmically interface with the web service.

***********