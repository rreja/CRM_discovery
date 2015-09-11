<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="meme.css">
    <![CDATA[
    /* START INCLUDED FILE "meme.css" */
        /* The following is the content of meme.css */
        body { background-color:white; font-size: 12px; font-family: Verdana, Arial, Helvetica, sans-serif;}

        div.help {
          display: inline-block;
          margin: 0px;
          padding: 0px;
          width: 12px;
          height: 13px;
          cursor: pointer;
          background-image: url(data:image/gif;base64,R0lGODlhDAANAIABANR0AP///yH5BAEAAAEALAAAAAAMAA0AAAIdhI8Xy22MIFgv1DttrrJ7mlGNNo4c+aFg6SQuUAAAOw==);
        }

        div.help:hover {
          background-image: url(data:image/gif;base64,R0lGODlhDAANAKEAANR0AP///9R0ANR0ACH+EUNyZWF0ZWQgd2l0aCBHSU1QACH5BAEAAAIALAAAAAAMAA0AAAIdDGynCe3PgoxONntvwqz2/z2K2ImjR0KhmSIZUgAAOw==);
        }
        
        p.spaced { line-height: 1.8em;}
        
        span.citation { font-family: "Book Antiqua", "Palatino Linotype", serif; color: #004a4d;}

        p.pad { padding-left: 30px; padding-top: 5px; padding-bottom: 10px;}

        td.jump { font-size: 13px; color: #ffffff; background-color: #00666a;
          font-family: Georgia, "Times New Roman", Times, serif;}

        a.jump { margin: 15px 0 0; font-style: normal; font-variant: small-caps;
          font-weight: bolder; font-family: Georgia, "Times New Roman", Times, serif;}

        h2.mainh {font-size: 1.5em; font-style: normal; margin: 15px 0 0;
          font-variant: small-caps; font-family: Georgia, "Times New Roman", Times, serif;}

        h2.line {border-bottom: 1px solid #CCCCCC; font-size: 1.5em; font-style: normal;
          margin: 15px 0 0; padding-bottom: 3px; font-variant: small-caps;
          font-family: Georgia, "Times New Roman", Times, serif;}

        h4 {border-bottom: 1px solid #CCCCCC; font-size: 1.2em; font-style: normal;
          margin: 10px 0 0; padding-bottom: 3px; font-family: Georgia, "Times New Roman", Times, serif;}

        h5 {margin: 0px}

        a.help { font-size: 9px; font-style: normal; text-transform: uppercase;
          font-family: Georgia, "Times New Roman", Times, serif;}

        div.pad { padding-left: 30px; padding-top: 5px; padding-bottom: 10px;}
        
        div.pad1 { margin: 10px 5px;}

        div.pad2 { margin: 25px 5px 5px;}
        h2.pad2 { padding: 25px 5px 5px;}

        div.pad3 { padding: 5px 0px 10px 30px;}

        div.box { border: 2px solid #CCCCCC; padding:10px; overflow: hidden;}

        div.bar { border-left: 7px solid #00666a; padding:5px; margin-top:25px; }

        div.subsection {margin:25px 0px;}

        img {border:0px none;}

        th.majorth {text-align:left;}
        th.minorth {font-weight:normal; text-align:left; width:8em; padding: 3px 0px;}
        th.actionth {font-weight:normal; text-align:left;}

        .block_td {height:25px;}
        .block_container {position:relative; width:98%; height:25px; padding:0px; margin: 0px 0px 0px 1em;}
        .block_motif {position:absolute; z-index:3; height:12px; top:0px; text-align:center; vertical-align:middle; background-color:cyan;}
        .block_rule {position:absolute; z-index:2; width:100%; height:1px; top:12px; left:0px; background-color:gray;}
        .block_plus_sym {position:absolute; z-index:4; line-height:12px; top:0px; left:-1em;}
        .block_minus_sym {position:absolute; z-index:4; line-height:12px; top:13px; left:-1em;}

        .tic_major {position:absolute; border-left:2px solid blue; height:0.5em; top:0em;}
        .tic_minor {position:absolute; border-left:1px solid blue; height:0.2em; top:0em;}
        .tic_label {position:absolute; top:0.5em;  height: 1em; text-align:center; vertical-align:middle}

        .explain h5 {font-size:1em; margin-left: 1em;}

        div.doc {margin-left: 2em; margin-bottom: 3em;}
        
        th.trainingset {
          border-bottom: thin dashed black; 
          font-weight:normal; 
          padding:0px 10px;
        }
        .dnaseq {
          font-weight: bold; 
          font-size: large; 
          font-family: 'Courier New', Courier, monospace;
        }
        .dna_A {
          color: rgb(204,0,0);
        }
        .dna_C {
          color: rgb(0,0,204);
        }
        .dna_G {
          color: rgb(255,179,0);
        }
        .dna_T {
          color: rgb(0,128,0);
        }
        div.pop_content {
          position:absolute;
          z-index:50;
          width:300px;
          padding: 5px;
          background: #E4ECEC;
          font-size: 12px;
          font-family: Arial;
          border-style: double;
          border-width: 3px;
          border-color: #AA2244;
          display:none;
        }
        div.pop_content > div.pop_close {
          float:right;
          bottom:0px;
        }
/*****************************************************************************
 * Program logo styling
 ****************************************************************************/
div.prog_logo {
  border-bottom: 0.25em solid #0f5f60;
  height: 4.5em;
  width: 24em;
  display:inline-block;
}
div.prog_logo img {
  float:left;
  width: 4em;
  border-style: none;
  margin-right: 0.2em;
}
div.prog_logo h1, div.prog_logo h1:hover, div.prog_logo h1:active, div.prog_logo h1:visited {
  margin:0;
  padding:0;
  font-family: Arial, Helvetica,  sans-serif;
  font-size: 3.2em;
  line-height: 1em;
  vertical-align: top;
  display: block;
  color: #026666;
  letter-spacing: -0.06em;
  text-shadow: 0.04em 0.06em 0.05em #666;
}
div.prog_logo h2, div.prog_logo h2:hover, div.prog_logo h2:active, div.prog_logo h2:visited {
  display: block;
  margin:0;
  padding:0;
  font-family: Helvetica, sans-serif;
  font-size: 0.9em;
  line-height: 1em;
  letter-spacing: -0.06em;
  color: black;
}

div.big.prog_logo {
  font-size: 18px;
}

    /* END INCLUDED FILE "meme.css" */
    ]]>
  </xsl:template>
</xsl:stylesheet>

