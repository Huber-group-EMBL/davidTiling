#!/usr/bin/perl -w

# perl cgi script to show gene transcript page from yeast tiling array analysis

use CGI;
use Symbol;
use AnyDBM_File ;
use Fcntl;

# introduce some variables:
my ($comment);
my @fields;
my ($name_col_no,$gene_col_no, $plot_col_no);
my ($start_col_no,$end_col_no, $strand_col_no);
my ($orf_col_no,$note_col_no,$chrom_col_no);
my ($chrom,$positionK,$positionUp, $positionDown);

# prepare head and tail of html-page
my $project_dir="/ebi/www/main/html/huber-srv/queryGene";

# directories where to find the samples' alongChromosome plots and definitions what these samples are:
my $baseDir = "http://www.ebi.ac.uk/huber-srv/queryGene/";
my @gene_pages_dir = ($baseDir."viz1", $baseDir."viz2", $baseDir."viz3", $baseDir."viz4", $baseDir."viz5");
my @sample_names   = ("Poly(A) Sample","Total-RNA Sample","Oligo(dT) reverse transcribed poly(A) sample (n=1)","Poly(A) 05_04_20","RNA-labeled Direct Hybridization (n=1)");
# the path to the alongChromosome plots to represent is put together out of these arrays

my $gffhashFileName = $project_dir."/"."gffhash"; # uncomment for web58-node1
#my $gffhashFileName = "gffhash"; # uncomment for local testing
my $cgiName = "http://www.ebi.ac.uk/huber-srv/cgi-bin/viewYeastTilingArray.pl";

# new version: load hash and get gene_pages from there
my %gffhash;

#dbmopen (%gffhash,$gffhashFileName,0222) || die "can't open hash: $!";
tie %gffhash, "AnyDBM_File", $gffhashFileName, O_RDONLY, 0
        or die "Cannot open file $gffhashFileName: $!\n";

$name_col_no = $gffhash{'Name'};
$gene_col_no = $gffhash{'gene'};
$plot_col_no = $gffhash{'plotfile'};
$start_col_no = $gffhash{'start'};
$end_col_no = $gffhash{'end'};
$strand_col_no = $gffhash{'strand'};
$orf_col_no = $gffhash{'orf_classification'};
$note_col_no = $gffhash{'Note'};
$chrom_col_no = $gffhash{'chr'};

#print " Column-Number: $name_col_no, $gene_col_no, $plot_col_no, $start_col_no, $end_col_no, $strand_col_no, $orf_col_no, $note_col_no, $chrom_col_no\n\n";
#my $numberOfGenePages = scalar keys %gffhash - 10;

# interpret gene argument:
my $q = CGI::new();


my $queryGene = $q->param("gene"); # uncomment for local testing
#my $queryGene = "YAR027W"; # comment for web58-node1
$queryGene =~ s/\s//g; # remove space characters in gene
$queryGene =~ tr/a-z/A-Z/; # convert to uppercase

my $querySamples="only2"; # set default here
# which samples should be shown? just Poly(A) and Total or all?
#if (defined($q->param("showSamples"))){
#  my $querySamples = $q->param("showSamples");
#}
 
if ($querySamples eq 'only2') {
  @gene_pages_dir = @gene_pages_dir[0,1];
  @sample_names = @sample_names[0,1];
}

# see if query gene exists:
my $queryStatus = 0; # initialize
my $queryOutput;

# Version 1.2: also allow chromosome coordinate:

# Part a: if the query begins with a character, it is possibly a gene name:
if (($queryGene =~ /^\D.*/) and ($queryGene =~ /^\w.*/)){

  # do we know that gene?
  if (defined($queryOutput =$gffhash{$queryGene})){

    $queryStatus = 1;
    @fields = split("\t",$queryOutput);

    # get chromosome position from filename, compute neighboring plots:
    $plotName = $fields[$plot_col_no];
    @nameComponents = split /_/,$plotName;
    $chrom = $nameComponents[0];
    $positionK = $nameComponents[1];
    if ($positionK > 0 ) { # there's nothing left of 0
	$positionDown =  sprintf "%04d", ($positionK - 5); # one picture left
    } else {
        $positionDown = $positionK;
    }
    $positionUp = sprintf "%04d", ($positionK + 5);  # one picture right

    # look up and insert according along-chromosome plots:  
    $queryOutput = "";
    for ($i=0; $i < @gene_pages_dir; $i++){
	$queryOutput = $queryOutput."<h2>".$sample_names[$i]."</h2>";
	# new version: set up arrows to link to adjacent picture
        # arrow left:
        $queryOutput = $queryOutput."<table><tr><td align=\"left\"><a href=\"".$cgiName."?gene=".$chrom.":".$positionDown."000\"><img src=\"".$baseDir."arrow_left.gif\" valign=\"middle\" border=\"0\"></a></td>";
        # picture itself:
	$queryOutput = $queryOutput."<td align=\"middle\">";
	$queryOutput = $queryOutput."<img src=\"".$gene_pages_dir[$i]."/".$fields[$plot_col_no].".jpg\" valign=\"middle\"></td>";
        # arrow right:
        $queryOutput = $queryOutput."<td align=\"right\"><a href=\"".$cgiName."?gene=".$chrom.":".$positionUp."000\"><img src=\"".$baseDir."arrow_right.gif\" valign=\"middle\" border=\"0\"></a></td></tr></table>";
    } # for each wanted sample

    $queryOutput = $queryOutput."<br /><center><h2>Annotation for $queryGene</h2></center>";
    $queryOutput = $queryOutput."<p><table><tr>\n"; # start header row
    $queryOutput = $queryOutput."<th>Gene Name</th>";
    $queryOutput = $queryOutput."<th>Alias(es)</th>";
    $queryOutput = $queryOutput."<th>Chromosome</th>";
    $queryOutput = $queryOutput."<th>Start [BP]</th>";
    $queryOutput = $queryOutput."<th>End [BP]</th>";
    $queryOutput = $queryOutput."<th>Strand</th>";
    $queryOutput = $queryOutput."<th>ORF-Classification</th>";
    $queryOutput = $queryOutput."<th>Notes</th>";
    $queryOutput = $queryOutput."</tr><tr align=\"center\">"; # start data row

    $queryOutput = $queryOutput."<td>".$fields[$name_col_no]."</td>"; 
    $queryOutput = $queryOutput."<td>".$fields[$gene_col_no]."</td>"; 
    $queryOutput = $queryOutput."<td>".$fields[$chrom_col_no]."</td>"; 
    $queryOutput = $queryOutput."<td>".$fields[$start_col_no]."</td>"; 
    $queryOutput = $queryOutput."<td>".$fields[$end_col_no]."</td>"; 
    $queryOutput = $queryOutput."<td>".$fields[$strand_col_no]."</td>"; 
    $queryOutput = $queryOutput."<td>".$fields[$orf_col_no]."</td>"; 
    $queryOutput = $queryOutput."<td>".$fields[$note_col_no]."</td>"; 
    $queryOutput = $queryOutput."</tr></table></p>\n"; 
  
  } else {  # query not defined as hash key

    $queryStatus = 0;

    # Version 1.1: partial matching
    my @mykeys = keys %gffhash;
    my @matches = grep(/$queryGene/i, @mykeys);
    my $number_of_matches = scalar(@matches);
    if ($number_of_matches < 1){
      $queryOutput = "Gene $queryGene not found! Please try again.";
    } elsif ($number_of_matches > 12){
      $queryOutput = "$number_of_matches gene names nearly matching $queryGene found! Please be more specific.";
    } else {
      $queryOutput = "Gene $queryGene not found! Did you mean one of these?<br>";
      $queryOutput = $queryOutput."<center><form action=\"$cgiName\">\n";
      $queryOutput = $queryOutput."<select name=\"gene\">";
      foreach $matchkey (@matches) {
	$queryOutput = $queryOutput."<option>$matchkey</option>";
      }
      $queryOutput = $queryOutput."</select>  ";
      $queryOutput = $queryOutput."<input type=\"submit\" value=\"Submit\">";
      $queryOutput = $queryOutput."</form></center>";
    }
  }

# Version 1.2: also allow chromosome coordinate
# b. possible chromosome coordinate:
} elsif ($queryGene =~ /^\d\d?:\d\d*$/){ # seems to be a chromosome coordinate?
  @fields = split(":",$queryGene);
  my $chrom = sprintf "%02d", $fields[0]; # format 3 to be shown as 03
  if ($chrom > 16){$queryOutput = "Yeast has only 16 chromosomes.";}
  else {
    use integer;
    # pictures computed in steps of 5 kB, so take lower 5kB step
    my $coord =  ($fields[1] / 5000) * 5;  # pictures computed in steps of 5 kB
    $coord = sprintf "%04d", $coord;
    # test if coordinate is still on chromosome;
    my @chrlengths = (230210, 813138, 316613, 1531914, 576869, 270148, 1090944, 562639, 439885, 745446, 666445, 1078173, 924430, 784328, 1091285, 948060, 85779); # taken from Bioconductors YEAST package
    $thischrlength = $chrlengths[($chrom-1)];
    if (($coord*1000) > $thischrlength){

      $queryOutput = "Chromosome $chrom is only $thischrlength bases long.";

    } else { # present according plots around this coordinate
      # get chromosome position from filename, compute neighboring plots:

	$positionK = $coord;
	if ($positionK > 0 ) { # there's nothing left of 0
	    $positionDown =  sprintf "%04d", ($positionK - 5); # one picture left
	} else {
	    $positionDown = $positionK;
	}
	$positionUp = sprintf "%04d", ($positionK + 5);  # one picture right

	$queryOutput = "";

	for ($i=0; $i<@gene_pages_dir; $i++){
	    $queryOutput = $queryOutput."<h2>$sample_names[$i]</h2>";
	    # new version: set up arrows to link to adjacent picture
	    # arrow left:
	    $queryOutput = $queryOutput."<table><tr><td align=\"left\"><a href=\"".$cgiName."?gene=".$chrom.":".$positionDown."000\"><img src=\"".$baseDir."arrow_left.gif\" valign=\"middle\" border=\"0\"></a></td>";
	    # picture itself:
	    $queryOutput = $queryOutput."<td align=\"middle\">";
	    $queryOutput = $queryOutput."<img src=\"".$gene_pages_dir[$i]."/".$chrom."_".$coord.".jpg\" valign=\"middle\"></td>";
	    # arrow right:
	    $queryOutput = $queryOutput."<td align=\"right\"><a href=\"".$cgiName."?gene=".$chrom.":".$positionUp."000\"><img src=\"".$baseDir."arrow_right.gif\" valign=\"middle\" border=\"0\"></a></td></tr></table>";
	    # old version:
	    #$queryOutput = $queryOutput."<p><img src=\"".$gene_pages_dir[$i]."/".$chrom."_".$coord.".jpg\"></p><br>";

       } # for each wanted sample

    } # else if plots do exist
  }# matches else { use integer;
} else {
  $queryOutput = "No result found for query $queryGene. Please try again!";
}

untie %gffhash;
#start html
print $q->header("text/html");

print <<EOF
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html lang="en">
<head>
<title>Tiling Array Transcriptome Analysis of Saccharomyces cerevisiae</title>
<meta http-equiv="Owner" content= "EMBL Outstation - Hinxton, European Bioinformatics Institute" />
<meta name="Author" content="Joern Toedling" />
<link rel="stylesheet" href="http://www.ebi.ac.uk/services/include/stylesheet.css" type="text/css" />
<link rel="shortcut icon" href="http://www.ebi.ac.uk/bookmark.ico" />
<meta http-equiv="content-type" content="text/html; charset=iso-8859-1" />
<script type="text/javascript"  src="http://www.ebi.ac.uk/include/master.js"></script>
<script type="text/javascript">
<!--
var emailaddress="toedling";
if (top != self){
    top.location.href = location.href;
   }
// -->
</script>
</head>
<body  marginwidth="0" marginheight="0" leftmargin="0" topmargin="0" rightmargin="0" bottommargin="0" onLoad="EbiPreloadImages('information');">

<!-- ############# see headers below for place to edit ################### -->

      <table width="100%" border="0" cellspacing="0" cellpadding="0"  class = "tabletop">
	  <tr>
	    <td width="270" height="65" align = "right"><a href="http://www.ebi.ac.uk/"><img  src="http://www.ebi.ac.uk/services/images/ebi_banner_1b.jpg" alt="EBI Home Page" width="270" height="65" border="0" /></a></td>
		<td valign = "top" align = "right" width = "100%"><table border = "0" cellspacing="0" cellpadding="0" class = "tabletop" width = "100%" height="65">
              <tr> 
                <td valign = "top" align = "right" colspan="2"><table border = "0" cellspacing="0" cellpadding="0" height = "28" class = "tablehead">
                    <tr> 
                      <td   class = "tablehead"  align = "left" valign = "bottom"><img alt="Image" src="http://www.ebi.ac.uk/services/images/top_corner.gif" width="28" height="28" /></td>
                        <form name = "Text1293FORM" action = "javascript:querySRS(document.forms[0].db[document.forms[0].db.selectedIndex].value, document.forms[0].qstr.value)" method = "post">
                        <td align = "center" valign = "middle"   class = "small" nowrap><span class = "smallwhite"><nobr>Get&nbsp;</nobr></span></td>
                        <td align = "center" valign = "middle"   class = "small"><span class = "small"><select  id = "FormsComboBox2" name = "db" class = "small">
                            <option value = "EMBL" selected >Nucleotide sequences</option>
							<option value = "SQUID">Protein sequences</option>
							<option value = "PDB">Protein structures</option>
							<option value = "INTERPRO">Protein signatures</option>
							<option value = "MEDLINE">Literature</option>
							<option value = "UNIPROT">Protein seq's [SRS]</option>
                          </select></span></td>
                        <td align = "center" valign = "middle"   class = "small" nowrap><span class = "smallwhite">&nbsp;for&nbsp;</span></td>
                        <td align = "center" valign = "middle"   class = "small"><span class = "small"><input id = "FormsEditField3" maxlength = "50" size = "7" name = "qstr"  class = "small" /></span></td>
                        <td align = "center" valign = "middle"   class = "small">&nbsp;</td>
                        <td align = "center" valign = "middle"   class = "small"><span class = "small"><input id = "FormsButton3" type = "submit" value = "Go" name = "FormsButton1" class = "small" /></span></td>
                        <td align = "center" valign = "middle"   class = "small" width="10" nowrap><a href = "#" class = "small2" onClick = "openWindow('http://www.ebi.ac.uk/help/DBhelp/dbhelp_frame.html'); return false;"><nobr>&nbsp;?&nbsp;</nobr></a></td>
                      </form>
                      <form name="google" action="http://www.google.com/u/ebi" method="get" onSubmit ="if (document.google.q.value=='') { alert('Please enter query.'); return false;}">
						<input type="hidden" name="hq" value="inurl:www.ebi.ac.uk" />
                        <td align = "center" valign = "middle"   class = "smallwhite" nowrap><span class = "smallwhite"><nobr>&nbsp;Site search&nbsp;</nobr></span></td>
                        <td align = "center" valign = "middle"   class = "small"><span class = "small"><input id="FormsEditField4" type="text" maxlength="50" size="7" name="q" class="small" /></span></td>
                        <td align = "center" valign = "middle"   class = "small">&nbsp;</td>
                        <td align = "center" valign = "middle"   class = "small"><span class = "small"><input id="FormsButton2" type="submit" value="Go" name="sa" class="small" /></span></td>
                        <td align = "center" valign = "middle"   class = "small" nowrap><nobr><a href = "#" class = "small2" onClick = "openWindow('http://www.ebi.ac.uk/help/help/sitehelp_frame.html'); return false;">&nbsp;?&nbsp;</a></nobr></td>
                      </form>
                    </tr>
                  </table></td>
              </tr>
              <tr> 
                <td align = "left" valign = "bottom"><a href="http://www.ebi.ac.uk/"><img  src="http://www.ebi.ac.uk/services/images/ebi_banner_2b.jpg" alt="EBI Home Page" width="66" height="29" border="0" /></a></td>
                <td align = "right" valign = "middle"><img alt="Image" src="http://www.ebi.ac.uk/services/images/thetopbar_b.gif" width="156" height="25" usemap="#Map" border="0" /></td>
              </tr>
            </table></td>
	</tr>
	<tr><td colspan = "2"><img alt="Image" src="http://www.ebi.ac.uk/services/images/trans.gif" width = "1" height = "5" /></td></tr>
</table> 
<table width="100%" border="0" cellspacing="0" cellpadding="0"  class = "tabletop" >
	<tr>
		<td width = "100%"><table width="679" border="0" cellspacing="0" cellpadding="0">
			<tr>
                <td width="97" height="18"><a href="http://www.ebi.ac.uk/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image8','','http://www.ebi.ac.uk/services/images/home_o.gif',1)"><img name="Image8" border="0" src="http://www.ebi.ac.uk/services/images/home.gif" width="97" height="18" /></a></td>
                <td width="97" height="18"><a href="http://www.ebi.ac.uk/Information/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image9','','http://www.ebi.ac.uk/services/images/about_o.gif',1)"><img name="Image9" border="0" src="http://www.ebi.ac.uk/services/images/about.gif" width="97" height="18" /></a></td>
                <td width="97" height="18"><a href="http://www.ebi.ac.uk/Groups/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image10','','http://www.ebi.ac.uk/services/images/research_o.gif',1)"><img name="Image10" border="0" src="http://www.ebi.ac.uk/services/images/research.gif" width="97" height="18" /></a></td>
				<td width="97" height="18"><a href="http://www.ebi.ac.uk/services/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image11','','http://www.ebi.ac.uk/services/images/services_o.gif',1)"><img name="Image11" border="0" src="http://www.ebi.ac.uk/services/images/services_o.gif" width="97" height="18" /></a></td>
                <td width="97" height="18"><a href="http://www.ebi.ac.uk/Tools/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image12','','http://www.ebi.ac.uk/services/images/utilities_o.gif',1)"><img name="Image12" border="0" src="http://www.ebi.ac.uk/services/images/utilities.gif" width="97" height="18" /></a></td>
                <td width="97" height="18"><a href="http://www.ebi.ac.uk/Databases/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image13','','http://www.ebi.ac.uk/services/images/databases_o.gif',1)"><img name="Image13" border="0" src="http://www.ebi.ac.uk/services/images/databases.gif" width="97" height="18" /></a></td>
                <td width="97" height="18"><a href="http://www.ebi.ac.uk/FTP/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image14','','http://www.ebi.ac.uk/services/images/downloads_o.gif',1)"><img name="Image14" border="0" src="http://www.ebi.ac.uk/services/images/downloads.gif" width="97" height="18" /></a></td> 
                <td width="97" height="18"><a href="http://www.ebi.ac.uk/Submissions/" onMouseOut="MM_swapImgRestore()" onMouseOver="MM_swapImage('Image15','','http://www.ebi.ac.uk/services/images/submissions_o.gif',1)"><img name="Image15" border="0" src="http://www.ebi.ac.uk/services/images/submissions.gif" width="97" height="18" /></a></td>
			</tr>
		</table></td></tr>
	<tr>
		<td width="100%" height = "5"  class = "tablehead" ><table width="100%" height = "5"  border="0" cellspacing="0" cellpadding="0">
              <tr> 
                <td width = "100%" height = "20" align = "center">
<!-- InstanceBeginEditable name="topnav" --><nobr><a href="http://www.ebi.ac.uk/huber-srv/queryGene/index.html" class = "white"><font size="+1">Tiling Array Transcriptome Analysis of Saccharomyces cerevisiae</font></a></nobr><!-- InstanceEndEditable --></td>
              </tr>
	      </table>
	    </td>
	  </tr>
      </table>
<!-- InstanceBeginEditable name="contents" -->
<!-- #################  Contents: edit here ######################## -->
      <center>
	<h2>New Query</h2>
        Enter gene symbol &nbsp;<i>(e.g.,  RPS16A, GIM3, GCN4)</i><br>
        OR chromosomal coordinate as chr:bp &nbsp;<i>(e.g., 13:552000)</i>
        <br><br>
 	<form action="$cgiName">
	  <input type="text" name="gene" SIZE=10 VALUE="$queryGene">
	  <input type="submit" value="Submit">
         </form>
        </center>
       <br>
       <hr>
EOF
;   
#          <br><br>
#          <input type="RADIO" name="showSamples" value="only2" CHECKED>
#          Poly(A)- and Total-RNA &nbsp;&nbsp;
#          <input type="RADIO" name="showSamples" value="all">
#          Include Additional Samples
#	</form>
#      </center>
#      <br>
#      <hr>
#EOF
#;
#}; # end of former multi-line comment

print $q->p("<center><h2>Transcription around $queryGene</h2></center>");

print $q->p("$queryOutput");

print <<EOF
<hr>
      <table width="100%">
	  <tr>
	    <td align="right">
	      <address><a href="mailto:toedling\@ebi.ac.uk">Joern Toedling</a></address>
	    </td>
	  </tr>
      </table>
EOF
;

print $q->end_html();



#
#   end of file
#
#
