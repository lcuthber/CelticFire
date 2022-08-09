#! /usr/bin/perl -w
use strict;
#author: Ulrike Löber (ulrike.loeber@mdc-berlin.de)
opendir my $dir_h, '.'
    or die "Cannot open directory: $!";

my @files = grep { /\.html$/ } readdir $dir_h;

closedir $dir_h;

open(OUTtax,">taxa_miga.txt") or die $!;
print OUTtax "ID\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";
open(OUTqual,">genome_quality_miga.txt") or die $!;
print OUTqual "ID\tcollection\tno_ess_genes\ttot_ess_genes\tcompleteness\tcontamination_perc\tquality\n";
open(OUTnovel,">tax_novelty_miga.txt") or die $!;
print OUTnovel "ID\tRoot\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tSubsepcies\tDataset\n";
open(OUTtaxpval,">tax_classification_pval_miga.txt") or die $!; # renamed before: tax_classification_miga.txt
print OUTtaxpval "ID\tdomain\tdpval\tphylum\tppval\tclass\tcpval\torder\topval\tfamily\tfpval\tgenus\tgpval\tspecies\tspval\n";
open(OUTanovel,">anovelty.r") or die $!; # renamed before: tax_classification_miga.txt
print OUTanovel "ID\tRoot\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tSubspecies\tDataset\n";


foreach my $file(@files){
my $numlin=0;
my $cont=0;
my $quality=0;
my $comp=0;

	open (INFILE, $file) or die $!;
	my @current=<INFILE>;
	close INFILE;

	foreach my $line(@current){
		chomp $line;
		$numlin++;
#<!DOCTYPE html>
#<html lang="en"><head>
#<meta http-equiv="content-type" content="text/html; charset=UTF-8"><style type="text/css">.turbolinks-progress-bar {
#  position: fixed;
#  display: block;
#  top: 0;
#  left: 0;
#  height: 3px;
#  background: #0076ff;
#  z-index: 9999;
#  transition: width 300ms ease-out, opacity 150ms 150ms ease-in;
#  transform: translate3d(0, 0, 0);
#}</style>
#  <!-- Global site tag (gtag.js) - Google Analytics -->
#  <script type="text/javascript" async="" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/analytics.js"></script><script async="" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/js.js"></script>
#  <script>
#    window.dataLayer = window.dataLayer || [];
#    function gtag(){dataLayer.push(arguments);}
#    gtag('js', new Date());
#    gtag('config', 'UA-43115844-2');
#  </script>
#  
#  <link rel="stylesheet" media="all" href="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/application-357e7d95314bd9c3e4e8601dfffd8c34bbf8a465f75720fb.css" data-turbolinks-track="true">
#  <script src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/application-29e8816e35782472f7087e4608067d071df540a952112b7d3.js" data-turbolinks-track="true"></script>
#  <script src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/plotly-latest.js" data-turbolinks-track="true" async="async"></script>
#  <!--[if lt IE 9]>
#   <script src="//cdnjs.cloudflare.com/ajax/libs/html5shiv/r29/html5.min.js">
#   </script>
#<![endif]-->
#
#  
#
#  
#<script charset="utf-8" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/momenttimelinetweet.js"></script><script charset="utf-8" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/timeline.js"></script><style id="plotly.js-style-global"></style><style id="plotly.js-style-modebar-e45a33"></style><style id="plotly.js-style-modebar-339cb0"></style><style id="plotly.js-style-modebar-3f488f"></style><style id="plotly.js-style-modebar-70ae2b"></style><style id="plotly.js-style-modebar-a5a950"></style><style id="plotly.js-style-modebar-eb0fda"></style><style id="plotly.js-style-modebar-96176f"></style><style id="plotly.js-style-modebar-f8c434"></style><style id="plotly.js-style-modebar-aa66f9"></style><style id="plotly.js-style-modebar-104ab7"></style><style id="plotly.js-style-modebar-2b5b1d"></style><style id="plotly.js-style-modebar-c282fe"></style><style id="plotly.js-style-modebar-1f4bb7"></style><style id="plotly.js-style-modebar-becd0a"></style><style id="plotly.js-style-modebar-564516"></style><style id="plotly.js-style-modebar-4da49c"></style><style id="plotly.js-style-modebar-ddf605"></style><style id="plotly.js-style-modebar-61e004"></style><style id="plotly.js-style-modebar-d7e934"></style><style id="plotly.js-style-modebar-b030dc"></style><style id="plotly.js-style-modebar-f331b5"></style><style id="plotly.js-style-modebar-330d8f"></style><style id="plotly.js-style-modebar-70e557"></style><style id="plotly.js-style-modebar-178ae0"></style><style id="plotly.js-style-modebar-21e13b"></style><style id="plotly.js-style-modebar-72b68e"></style><style id="plotly.js-style-modebar-6316ac"></style><style id="plotly.js-style-modebar-c8c987"></style><style id="plotly.js-style-modebar-adfd0a"></style><style id="plotly.js-style-modebar-f5c06a"></style><style id="plotly.js-style-modebar-49fa71"></style><style id="plotly.js-style-modebar-6244f0"></style><style id="plotly.js-style-modebar-146764"></style><style id="plotly.js-style-modebar-8112f0"></style><style id="plotly.js-style-modebar-232619"></style><style id="plotly.js-style-modebar-898e2e"></style><style id="plotly.js-style-modebar-c690a4"></style><style id="plotly.js-style-modebar-e810c1"></style><style id="plotly.js-style-modebar-0325e4"></style><style id="plotly.js-style-modebar-04879f"></style><style id="plotly.js-style-modebar-47c2a8"></style><style id="plotly.js-style-modebar-7685ab"></style><style id="plotly.js-style-modebar-bcfe29"></style><style id="plotly.js-style-modebar-d1be2a"></style><style id="plotly.js-style-modebar-73a0f0"></style><style id="plotly.js-style-modebar-b1343b"></style><style id="plotly.js-style-modebar-f47e41"></style><style id="plotly.js-style-modebar-ba2719"></style><style id="plotly.js-style-modebar-715c3f"></style><style id="plotly.js-style-modebar-01c94c"></style><style id="plotly.js-style-modebar-3f1bbd"></style><style id="plotly.js-style-modebar-b11ce6"></style><style id="plotly.js-style-modebar-c7e799"></style><style id="plotly.js-style-modebar-d076dc"></style><style id="plotly.js-style-modebar-cb080f"></style><style id="plotly.js-style-modebar-da1ea6"></style><style id="plotly.js-style-modebar-399080"></style><style id="plotly.js-style-modebar-d68bd6"></style><style id="plotly.js-style-modebar-ee4193"></style><style id="plotly.js-style-modebar-27bc92"></style><style id="plotly.js-style-modebar-e811c5"></style><style id="plotly.js-style-modebar-3f7dee"></style><style id="plotly.js-style-modebar-77057e"></style><style id="plotly.js-style-modebar-87b1dd"></style><style id="plotly.js-style-modebar-61b544"></style><style id="plotly.js-style-modebar-5ede3f"></style><style id="plotly.js-style-modebar-325a61"></style><style id="plotly.js-style-modebar-2425c1"></style><style id="plotly.js-style-modebar-0619c3"></style><style id="plotly.js-style-modebar-fbf014"></style><style id="plotly.js-style-modebar-275d94"></style><style id="plotly.js-style-modebar-eea733"></style><style id="plotly.js-style-modebar-8beb01"></style><style id="plotly.js-style-modebar-c2f1b3"></style><style id="plotly.js-style-modebar-f1f27d"></style><style id="plotly.js-style-modebar-7fd8f0"></style><style id="plotly.js-style-modebar-57a10f"></style><style id="plotly.js-style-modebar-c2631e"></style><style id="plotly.js-style-modebar-baa47b"></style><style id="plotly.js-style-modebar-a02551"></style><style id="plotly.js-style-modebar-786c6e"></style><style id="plotly.js-style-modebar-6dc0ab"></style><style id="plotly.js-style-modebar-143918"></style><style id="plotly.js-style-modebar-5ece87"></style><style id="plotly.js-style-modebar-db23a5"></style><style id="plotly.js-style-modebar-1f55e0"></style><style id="plotly.js-style-modebar-d9ae03"></style><style id="plotly.js-style-modebar-b2053e"></style><style id="plotly.js-style-modebar-df7a5a"></style><style id="plotly.js-style-modebar-8bf845"></style><style id="plotly.js-style-modebar-06c041"></style><style id="plotly.js-style-modebar-7820de"></style><style id="plotly.js-style-modebar-2a9519"></style><style id="plotly.js-style-modebar-452499"></style><style id="plotly.js-style-modebar-701e14"></style><style id="plotly.js-style-modebar-a9d988"></style><style id="plotly.js-style-modebar-c6a40a"></style><style id="plotly.js-style-modebar-f1997e"></style><style id="plotly.js-style-modebar-96fbb5"></style><style id="plotly.js-style-modebar-a17b43"></style><style id="plotly.js-style-modebar-ab59e0"></style><style id="plotly.js-style-modebar-7c69f0"></style><style id="plotly.js-style-modebar-0918bc"></style><style id="plotly.js-style-modebar-d83999"></style><style id="plotly.js-style-modebar-969589"></style><style id="plotly.js-style-modebar-7fc144"></style><style id="plotly.js-style-modebar-1fdd31"></style><style id="plotly.js-style-modebar-37e71e"></style><style id="plotly.js-style-modebar-da299f"></style><style id="plotly.js-style-modebar-bbcf3a"></style><style id="plotly.js-style-modebar-a0a63f"></style><style id="plotly.js-style-modebar-d254ad"></style><style id="plotly.js-style-modebar-30c761"></style><style id="plotly.js-style-modebar-65086d"></style><style id="plotly.js-style-modebar-98945d"></style><style id="plotly.js-style-modebar-7de5e6"></style><style id="plotly.js-style-modebar-fabe0b"></style><style id="plotly.js-style-modebar-099047"></style><style id="plotly.js-style-modebar-a5f986"></style><style id="plotly.js-style-modebar-add223"></style><style id="plotly.js-style-modebar-aa6e10"></style><style id="plotly.js-style-modebar-b42d6b"></style><style id="plotly.js-style-modebar-9ccd77"></style><style id="plotly.js-style-modebar-e7cdf0"></style><style id="plotly.js-style-modebar-0875ca"></style><style id="plotly.js-style-modebar-710eed"></style><style id="plotly.js-style-modebar-571de4"></style><style id="plotly.js-style-modebar-11ba89"></style><style id="plotly.js-style-modebar-090bb3"></style><style id="plotly.js-style-modebar-eaaa19"></style><style id="plotly.js-style-modebar-610a48"></style><style id="plotly.js-style-modebar-9373b2"></style><style id="plotly.js-style-modebar-aa5895"></style><style id="plotly.js-style-modebar-e3214f"></style><style id="plotly.js-style-modebar-61a99b"></style><style id="plotly.js-style-modebar-49ba5b"></style><style id="plotly.js-style-modebar-ff0381"></style><style id="plotly.js-style-modebar-f113f3"></style><style id="plotly.js-style-modebar-a65f83"></style><style id="plotly.js-style-modebar-a713f6"></style><style id="plotly.js-style-modebar-c8b647"></style><style id="plotly.js-style-modebar-e477b6"></style><style id="plotly.js-style-modebar-51306f"></style><style id="plotly.js-style-modebar-8bfe64"></style><style id="plotly.js-style-modebar-1b47d7"></style><style id="plotly.js-style-modebar-dfcb33"></style><style id="plotly.js-style-modebar-5e5de0"></style><style id="plotly.js-style-modebar-fdd5ad"></style><style id="plotly.js-style-modebar-7381c2"></style><style id="plotly.js-style-modebar-1b6fc1"></style><style id="plotly.js-style-modebar-092c82"></style><style id="plotly.js-style-modebar-365a10"></style><style id="plotly.js-style-modebar-a53d25"></style><style id="plotly.js-style-modebar-142656"></style><style id="plotly.js-style-modebar-0ba098"></style><style id="plotly.js-style-modebar-7abb57"></style><style id="plotly.js-style-modebar-51f7db"></style><style id="plotly.js-style-modebar-7f5ed6"></style><style id="plotly.js-style-modebar-576646"></style><style id="plotly.js-style-modebar-4284fc"></style><style id="plotly.js-style-modebar-8181eb"></style><style id="plotly.js-style-modebar-664dd9"></style><style id="plotly.js-style-modebar-73c59a"></style><style id="plotly.js-style-modebar-d895e7"></style><style id="plotly.js-style-modebar-e7cfd0"></style><style id="plotly.js-style-modebar-252143"></style><style id="plotly.js-style-modebar-67eaae"></style><style id="plotly.js-style-modebar-7b4993"></style><style id="plotly.js-style-modebar-535af3"></style><style id="plotly.js-style-modebar-847d5b"></style><style id="plotly.js-style-modebar-d09e92"></style><style id="plotly.js-style-modebar-b34ca9"></style><style id="plotly.js-style-modebar-e1b91e"></style><style id="plotly.js-style-modebar-6e7beb"></style><style id="plotly.js-style-modebar-c88d8d"></style><style id="plotly.js-style-modebar-7be723"></style><style id="plotly.js-style-modebar-285e28"></style><style id="plotly.js-style-modebar-ead4f2"></style><title>Streptococcus mitis 27098 8 93 | MiGA Online</title><meta name="csrf-param" content="authenticity_token"><meta name="csrf-token" content="ittOTgNXN1ooYj0cKNHVCaz/cnK86N3+vRpA7xUnSFUx+SVMcQvkJwXGlK3L95EOHy5n78EoyFEcdA/wilxnvg=="><style id="plotly.js-style-modebar-6bf225"></style></head>
#<body>
#  <header class="navbar navbar-fixed-top navbar-inverse">
#  <div class="container">
#    <nav>
#      <a id="logo" href="http://microbial-genomes.org/">MiGA</a>
#      <div class="crumb">
#          
#          / <a href="http://microbial-genomes.org/projects">Projects</a>
#          / <a href="http://microbial-genomes.org/projects/20">TypeMat</a>
#          
#          
#          
#          
#          
#          
#          
#        / Streptococcu...
#      </div>
#      <ul class="nav navbar-nav navbar-right">
#        <li>
#          <a class="dropdown-toggle" data-toggle="dropdown" type="button">
#            Projects
#            <span class="caret"></span>
#          </a>
#          <ul class="dropdown-menu" role="menu" aria-labelledby="dropdownMenu">
#              <li>
#                <a href="http://microbial-genomes.org/projects">
#                  <i class="glyphicon glyphicon-globe"> </i> Public
#</a>              </li>
#              <li>
#                <a href="http://microbial-genomes.org/projects?user_contributed=true">
#                  <i class="glyphicon glyphicon-user"> </i> User-contributed
#</a>              </li>
#              <li>
#                <a href="http://microbial-genomes.org/projects?type=clade">
#                  <i class="glyphicon glyphicon-tree-deciduous"> </i> Clade
#</a>              </li>
#          </ul>
#        </li>
#            <li>
#              <a href="http://microbial-genomes.org/query_datasets?complete_new=true">
#                <span class="badge">88 ready</span>
#</a>            </li>
#          <li>
#            <a class="dropdown-toggle" data-toggle="dropdown" type="button">
#              <img alt="Ulrike Löber" class="gravatar img-circle" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/556098948d99e1dcacfd718a6fe0f58b.png">
#              <span class="caret"></span>
#            </a>
#            <ul class="dropdown-menu" role="menu" aria-labelledby="dropdownMenu">
#              <li><a href="http://microbial-genomes.org/dashboard">
#                <i class="glyphicon glyphicon-dashboard"> </i> Dashboard
#              </a></li>
#              <li><a href="http://microbial-genomes.org/query_datasets">
#                <i class="glyphicon glyphicon-list"> </i> Query datasets
#              </a></li>
#              <li class="divider"></li>
#              <li><a rel="nofollow" data-method="delete" href="http://microbial-genomes.org/logout">
#                <i class="glyphicon glyphicon-user"> </i> Log out
#              </a></li>
#            </ul>
#          </li>
#      </ul>
#    </nav>
#  </div>
#</header>
#
#  <div class="container wrapper">
#    <div class="row">
#  <aside class="col-md-4">
#    <section class="dataset-title">
#      <div class="dataset-item query-dataset">
#  <a class="label label-primary pull-left project-code" style="background:hsl(232.82, 55.48%, 55.74%);" title="TypeMat" href="http://microbial-genomes.org/projects/20">TM</a>
#  <span class="dataset-name">
#        <i class="glyphicon glyphicon-ok-circle small text-light-success" rel="tooltip" title="" data-original-title="ready"> </i>
#        <i class="glyphicon glyphicon-adjust small text-light-success" rel="tooltip" title="" data-original-title="Excellent quality genome"> </i>
#    <b><a href="http://microbial-genomes.org/query_datasets/M:-WH_O5U">Streptococcus mitis 27098 8 93</a></b>
if ($line=~m/<b><a href=\"http:\/\/microbial-genomes.org\/query_datasets\//){
	$line=~s/.*\">//g;
	$line=~s/<.*//g;
	$line=~s/ /_/g;
	print OUTtax "$line\t";
	print OUTqual "$line\t";
	print OUTnovel "$line\t";
	print OUTtaxpval "$line\t";
}

#  </span>
#  <span class="timestamp">
#    Created about 1 month ago
#  </span>
#</div>
#
#
#    </section>
#    <section class="dataset-item">
#      <div class="dataset-item reference-dataset">
#  <div class="panel-group" id="dataset-metadata" role="tablist">
#    <div class="panel panel-default"><div class="panel-heading" role="tab"><h4 class="panel-title"><a class="btn" data-toggle="collapse" data-parent="#dataset-metadata" href="#dataset-metadata-taxonomy">Taxonomy</a></h4></div><div class="panel-collapse collapse in" id="dataset-metadata-taxonomy"><div class="panel-body">
#      <ul class="panel-body taxonomy">
#        <li>
#          <span class="text-muted">d:</span>
#            <i>Bacteria</i>
#        </li>
#        <li>
#          <span class="text-muted">p:</span>
#            <i>Firmicutes</i>
#        </li>
#        <li>
#          <span class="text-muted">c:</span>
#            <i>Bacilli</i>
#        </li>
#        <li>
#          <span class="text-muted">o:</span>
#            <i>Lactobacillales</i>
#        </li>
#        <li>
#          <span class="text-muted">f:</span>
#            <i>Streptococcaceae</i>
#        </li>
#        <li>
#          <span class="text-muted">g:</span>
#            <i>Streptococcus</i>
#        </li>
#      </ul>
#</div></div></div>
#
#    <div class="panel panel-default"><div class="panel-heading" role="tab"><h4 class="panel-title"><a class="btn" data-toggle="collapse" data-parent="#dataset-metadata" href="#dataset-metadata-type">Type</a></h4></div><div class="panel-collapse collapse " id="dataset-metadata-type"><div class="panel-body">
#      <b>Genome</b>: The genome from an isolate
#</div></div></div>  
#</div>
#
#</div>
#
#
#    </section>
#    <section>
#      <h4>Settings</h4>
#      <a class="btn btn-block btn-default" rel="nofollow" data-method="post" href="http://microbial-genomes.org/query_datasets/M:-WH_O5U/unread">Mark dataset as unseen</a>
#      <a class="btn btn-block btn-warning" data-confirm="Are you sure? This action cannot be undone." rel="nofollow" data-method="delete" href="http://microbial-genomes.org/query_datasets/M:-WH_O5U">Remove dataset</a>
#    </section>
#  </aside>
#  <div class="col-md-8">
#
#
#
#        <div id="query_dataset_result_stats_42005"></div>
#        <div id="query_dataset_result_distances_42005">
#
#<h3 role="button" data-toggle="collapse" class="result-title" id="result-distances" href="#result-cnt-distances" aria-expanded="true" aria-controls="result-cnt-distances">
#  <i class="glyphicon glyphicon-tasks"> </i>
#  Distances
#</h3>
#
#<div class="result collapse in" id="result-cnt-distances" aria-expanded="true">
#<div class="result-body">
#    <p>
#    The closest relatives found by MiGA in the database were
#    <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_NCTC_12261_NZ_CP028414"><span style="display:inline;">Streptococcus mitis NCTC 12261 NZ CP028414</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#    (95.46% ANI) and
#    <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_GCA_002014755"><span style="display:inline;">Streptococcus mitis GCA 002014755</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#    (95.43% ANI).
#</p>
#  <h4>Taxonomic classification
#  <sup class="info-msg-button"><i class="glyphicon glyphicon-info-sign text-info" data-toggle="modal" data-target="#info-msg-96387d7d-4377-4300-ac3b-b7014a337338"> </i></sup></h4>
#  <p>The dataset most likely belongs to the <b>genus</b> <i class="tax-name">Streptococcus</i> (p-value: 0.0047) and probably belongs to the <b>species</b> <i class="tax-name">Streptococcus mitis</i> (p-value: 0.08). </p><div class="small comment"><div class="taxonomy-tree"><b><span class="badge">root</span>  (p-value 0****)</b><div class="taxonomy-tree"><b><span class="badge">domain</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Bacteria%22">Bacteria</a></i> (p-value 0****)</b><div class="taxonomy-tree"><b><span class="badge">phylum</span> <span class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Firmicutes%22">Firmicutes</a></span> (p-value 0****)</b><div class="taxonomy-tree"><b><span class="badge">class</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Bacilli%22">Bacilli</a></i> (p-value 0****)</b><div class="taxonomy-tree"><b><span class="badge">order</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Lactobacillales%22">Lactobacillales</a></i> (p-value 0.000216****)</b><div class="taxonomy-tree"><b><span class="badge">family</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcaceae%22">Streptococcaceae</a></i> (p-value 0.000756****)</b><div class="taxonomy-tree"><b><span class="badge">genus</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcus%22">Streptococcus</a></i> (p-value 0.0047****)</b><div class="taxonomy-tree"><b><span class="badge">species</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcus+mitis%22">Streptococcus mitis</a></i> (p-value 0.0798**)</b><div class="taxonomy-tree"><b class="text-muted"><span class="badge">subspecies</span>  (p-value 0.703)</b><div class="taxonomy-tree"><b class="text-muted"><span class="badge">dataset</span> <span class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcus+mitis+NCTC+12261%22">Streptococcus mitis NCTC 12261</a></span> (p-value 0.738)</b></div></div></div></div></div></div></div></div></div></div><br><span class="text-muted">Significance at p-value below: *0.5, **0.1, ***0.05, ****0.01</span></div><p></p>
	elsif($line=~m/\<div class=\"taxonomy-tree\"\>\<b\>\<span class=\"badge\"\>/){
#		print "$line\n";
		my @array=split ("</b><div",$line);
		my $rootpval=$array[0];
		$rootpval=~s/.*p-value //;
		$rootpval=~s/\*.*//;
		$rootpval=~s/\).*//;
		print OUTanovel "$rootpval\t";
##<div class="taxonomy-tree"><b><span class="badge">domain</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Bacteria%22">Bacteria</a></i> (p-value 0****)</b>
		my $domain=$array[1];
		$domain=~s/<\/a>.*//;
		$domain=~s/.*>//;
		my $pval=$array[1];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTtax "$domain\t";
		print OUTtaxpval "$domain\t$pval\t";
		print OUTanovel "$pval\t";
##<div class="taxonomy-tree"><b><span class="badge">phylum</span> <span class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Firmicutes%22">Firmicutes</a></span> (p-value 0****)</b>
		my $phylum=$array[2];
		$phylum=~s/<\/a>.*//;
		$phylum=~s/.*>//;
		$pval=$array[2];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTtax "$phylum\t";
		print OUTtaxpval "$phylum\t$pval\t";
		print OUTanovel "$pval\t";
##<div class="taxonomy-tree"><b><span class="badge">class</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Bacilli%22">Bacilli</a></i> (p-value 0****)</b>
		my $class=$array[3];
		$class=~s/<\/a>.*//;
		$class=~s/.*>//;
		$pval=$array[3];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTtax "$class\t";
		print OUTtaxpval "$class\t$pval\t";
		print OUTanovel "$pval\t";

##<div class="taxonomy-tree"><b><span class="badge">order</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Lactobacillales%22">Lactobacillales</a></i> (p-value 0.000216****)</b>
		my $order=$array[4];
		$order=~s/<\/a>.*//;
		$order=~s/.*>//;
		$pval=$array[4];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTtax "$order\t";
		print OUTtaxpval "$order\t$pval\t";
		print OUTanovel "$pval\t";

##<div class="taxonomy-tree"><b><span class="badge">family</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcaceae%22">Streptococcaceae</a></i> (p-value 0.000756****)</b>
		my $family=$array[5];
		$family=~s/<\/a>.*//;
		$family=~s/.*>//;
		$pval=$array[5];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		if($array[6]=~m/Gemella/){	#Gemella has no family assigned -> missing in MIGA annotation (checked ncbi tax browser 06.02.2021)
			$family="NA";
		}
		print OUTtax "$family\t";
		print OUTtaxpval "$family\t$pval\t";
		print OUTanovel "$pval\t";

##<div class="taxonomy-tree"><b><span class="badge">genus</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcus%22">Streptococcus</a></i> (p-value 0.0047****)</b>
		my $genus=$array[6];
		$genus=~s/<\/a>.*//;
		$genus=~s/.*>//;
		$pval=$array[6];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTtax "$genus\t";
		print OUTtaxpval "$genus\t$pval\t";
		print OUTanovel "$pval\t";

##<div class="taxonomy-tree"><b><span class="badge">species</span> <i class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcus+mitis%22">Streptococcus mitis</a></i> (p-value 0.0798**)</b>
		my $species=$array[7];
		$species=~s/<\/a>.*//;
		$species=~s/.*>//;
		$pval=$array[7];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTtax "$species";
		print OUTtaxpval "$species\t$pval\t";
		print OUTanovel "$pval\t";

##<div class="taxonomy-tree"><b class="text-muted"><span class="badge">subspecies</span>  (p-value 0.703)</b>
		my $subspecies=$array[8];
		$subspecies=~s/<\/a>.*//;
		$subspecies=~s/.*>//;
		$pval=$array[8];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTtaxpval"$subspecies\t$pval";

		print OUTtax "\n";
		print OUTtaxpval "\n";
		print OUTanovel "$pval\t";
##<div class="taxonomy-tree"><b class="text-muted"><span class="badge">dataset</span> <span class="tax-name"><a href="http://microbial-genomes.org/projects/20/search?q=tax%3A%22Streptococcus+mitis+NCTC+12261%22">Streptococcus mitis NCTC 12261</a></span> (p-value 0.738)</b></div></div></div></div></div></div></div></div></div></div><br><span class="text-muted">Significance at p-value below: *0.5, **0.1, ***0.05, ****0.01</span></div><p></p>
		my $dataset=$array[9];
		$pval=$array[9];
		$pval=~s/.*p-value //;
		$pval=~s/\*.*//;
		$pval=~s/\).*//;
		print OUTanovel "$pval\n";
	}
#  <h4>Taxonomic novelty
#  <sup class="info-msg-button"><i class="glyphicon glyphicon-info-sign text-info" data-toggle="modal" data-target="#info-msg-fa40afad-5f94-4847-aa94-f7ecbb85b004"> </i></sup></h4>
#  <p>The dataset most likely belongs to a <b>subspecies</b> not represented in the database (p-value: 0.00055), highest taxonomic rank with p-value ≤ 0.01. It probably belongs to a <b>species</b> not represented in the database (p-value: 0.026), highest taxonomic rank with p-value ≤ 0.1. </p><div class="text-muted small comment"><b>P-values:</b> <b>root</b> 1, <b>domain</b> 1, <b>phylum</b> 0.999, <b>class</b> 0.998, <b>order</b> 0.993, <b>family</b> 0.983, <b>genus</b> 0.954, <b>species</b> 0.0256, <b>subspecies</b> 0.000545, <b>dataset</b> 0.000206.</div><p></p>
	if ($line=~m/<p>The dataset most likely belongs to a/){
		my @temp=split("<b>",$line);
##<b>root</b> 1, 
		my $root=$temp[4];
		$root=~s/.*<b>//;
		my $val=$root;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>domain</b> 1, 
		my $domain=$temp[5];
		$domain=~s/.*<b>//;
		$val=$domain;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>phylum</b> 0.999, 
		my $phylum=$temp[6];
		$phylum=~s/.*<b>//;
		$val=$phylum;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>class</b> 0.998, 
		my $class=$temp[7];
		$class=~s/.*<b>//;
		$val=$class;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>order</b> 0.993, 
		my $order=$temp[8];
		$order=~s/.*<b>//;
		$val=$order;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>family</b> 0.983, 
		my $family=$temp[9];
		$family=~s/.*<b>//;
		$val=$family;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>genus</b> 0.954, 
		my $genus=$temp[10];
		$genus=~s/.*<b>//;
		$val=$genus;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>species</b> 0.0256, 
		my $species=$temp[11];
		$species=~s/.*<b>//;
		$val=$species;
		$val=~s/.*> //;
		$val=~s/,.*//;
		print OUTnovel "$val\t";
##<b>subspecies</b> 0.000545, 
		my $subspecies=$temp[12];
		$subspecies=~s/.*<b>//;
		$val=$subspecies;
		$val=~s/.*> //;
		$val=~s/,.*//;
		if  ($subspecies=~m/<\/div/){
			$val=~s/\.<.*//;		
			print OUTnovel "$val\tNA\n";
		}
		else{
			print OUTnovel "$val\t";
			my $dataset=$temp[13];
			$dataset=~s/.*<b>//;
			$val=$dataset;
			$val=~s/.*> //;
			$val=~s/\.<.*//;
			print OUTnovel "$val\n";
		}
##<b>dataset</b> 0.000206.</div><p></p>
	}
#<h4>Genome relatedness
#<sup class="info-msg-button"><i class="glyphicon glyphicon-info-sign text-info" data-toggle="modal" data-target="#info-msg-0d90c62d-ffa8-40ce-9b18-a9609fbd64d4"> </i></sup></h4>
#  <p>
#    This dataset was placed in the context of some representative genomes from
#    the reference collection. See below the <b>ref tree</b> (in Newick format)
#    and <b>ref tree pdf</b> (in PDF).
#  </p>
#<div class="panel-group" id="distances" role="tablist" aria-multiselectable="true">
#  <div class="panel panel-default">
#    <div class="panel-heading" role="tab" id="ani-h">
#      <h4 class="panel-title">
#        <a role="button" data-toggle="collapse" data-parent="#distances" class="btn" href="#ani-b" aria-expanded="false" aria-controls="ani-b">
#          ANI table
#        </a>
#      </h4>
#    </div>
#    <div id="ani-b" class="panel-collapse collapse" role="tabpanel" aria-labelledby="ani-h">
#      <div class="panel-body table-responsive">
#        <table class="table table-hover">
#          <thead><tr>
#            <th>Dataset</th>
#            <th>ANI (%)</th>
#            <th>Standard deviation (ANI%)</th>
#            <th>Fraction of genome
#              shared (%)</th>
#          </tr></thead>
#          <tbody>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_NCTC_12261_NZ_CP028414"><span style="display:inline;">Streptococcus mitis NCTC 12261 NZ CP028414</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="95.4569">95.46</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="86.49517684887459">86.5</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_GCA_002014755"><span style="display:inline;">Streptococcus mitis GCA 002014755</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="95.4344">95.43</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="87.21311475409836">87.21</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_GCA_900459425"><span style="display:inline;">Streptococcus mitis GCA 900459425</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="95.3942">95.39</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="87.68">87.68</span>
#                </td>
#              </tr>
#          </tbody>
#        </table>
#      </div>
#    </div>
#  </div>
#  <div class="panel panel-default">
#    <div class="panel-heading" role="tab" id="aai-h">
#      <h4 class="panel-title">
#        <a role="button" data-toggle="collapse" data-parent="#distances" class="btn" href="#aai-b" aria-expanded="false" aria-controls="aai-b">
#          AAI table
#        </a>
#      </h4>
#    </div>
#    <div id="aai-b" class="panel-collapse collapse" role="tabpanel" aria-labelledby="aai-h">
#      <div class="panel-body table-responsive">
#        <table class="table table-hover">
#          <thead><tr>
#            <th>Dataset</th>
#            <th>AAI (%)</th>
#            <th>Standard deviation (AAI%)</th>
#            <th>Fraction of proteins
#              shared (%)</th>
#          </tr></thead>
#          <tbody>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_GCA_900459425"><span style="display:inline;">Streptococcus mitis GCA 900459425</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="96.61707632600272">96.62</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="7.268760763934852">7.27</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="89.36416184971098">89.36</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_NCTC_12261_NZ_CP028414"><span style="display:inline;">Streptococcus mitis NCTC 12261 NZ CP028414</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="96.59805825242724">96.6</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="7.268793704969919">7.27</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="88.99769585253456">89.0</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mitis_GCA_002014755"><span style="display:inline;">Streptococcus mitis GCA 002014755</span><sup title="assembly designated as neotype" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="96.59076227390211">96.59</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="7.291304770662442">7.29</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="89.63520555877244">89.64</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_chosunense_GCA_003626515"><span style="display:inline;">Streptococcus chosunense GCA 003626515</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="95.0616456536165">95.06</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="9.653183772556982">9.65</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="86.16352201257861">86.16</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_gwangjuense_NZ_CP032621"><span style="display:inline;">Streptococcus gwangjuense NZ CP032621</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="94.81453947368418">94.81</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="10.49347819630664">10.49</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="84.11732152739347">84.12</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_pseudopneumoniae_ATCC_BAA_960___CCUG_49455_GCA_002087075"><span style="display:inline;">Streptococcus pseudopneumoniae ATCC BAA 960   CCUG 49455 GCA 002087075</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="93.67004701141697">93.67</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="10.229012456436463">10.23</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="82.40177089097952">82.4</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_pneumoniae_GCA_001679535"><span style="display:inline;">Streptococcus pneumoniae GCA 001679535</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="93.287266435986">93.29</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="10.11203968084144">10.11</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="79.96679579413393">79.97</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_oralis_ATCC_35037_GCA_000148565"><span style="display:inline;">Streptococcus oralis ATCC 35037 GCA 000148565</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="89.83377808988762">89.83</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="12.202387558610681">12.2</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="78.80464858882125">78.8</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_oralis_ATCC_35037_GCA_000164095"><span style="display:inline;">Streptococcus oralis ATCC 35037 GCA 000164095</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="89.73352033660589">89.73</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="12.427443778004967">12.43</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="79.266259032796">79.27</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_infantis_ATCC_700779_GCA_000260755"><span style="display:inline;">Streptococcus infantis ATCC 700779 GCA 000260755</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="86.3796095444684">86.38</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="13.081916807294489">13.08</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="76.53569452130603">76.54</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_cristatus_ATCC_51100_GCA_000222765"><span style="display:inline;">Streptococcus cristatus ATCC 51100 GCA 000222765</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="75.68498475609753">75.68</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="16.711476317219507">16.71</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="72.606530160487">72.61</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_gordonii_NZ_LS483341"><span style="display:inline;">Streptococcus gordonii NZ LS483341</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="74.67886792452843">74.68</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="17.17614950813863">17.18</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="73.32595462091865">73.33</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_australis_ATCC_700641_GCA_000222745"><span style="display:inline;">Streptococcus australis ATCC 700641 GCA 000222745</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="74.4355855855857">74.44</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="17.571124732263876">17.57</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="73.71333702268954">73.71</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_sanguinis_NZ_LS483385"><span style="display:inline;">Streptococcus sanguinis NZ LS483385</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="74.2908951798011">74.29</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="16.88324916179028">16.88</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="72.32982844493635">72.33</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_anginosus_NZ_LR134283"><span style="display:inline;">Streptococcus anginosus NZ LR134283</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="74.11491935483875">74.11</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="16.05314344636658">16.05</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="68.62202545655784">68.62</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_parasanguinis_GCA_900459355"><span style="display:inline;">Streptococcus parasanguinis GCA 900459355</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="73.99746457867269">74.0</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="17.258830405897747">17.26</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="74.21140011068069">74.21</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_intermedius_NZ_LS483436"><span style="display:inline;">Streptococcus intermedius NZ LS483436</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="73.47399836467694">73.47</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="16.484338640025495">16.48</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="67.68123962368567">67.68</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_himalayensis_NZ_CP016953"><span style="display:inline;">Streptococcus himalayensis NZ CP016953</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="71.76789215686289">71.77</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="17.500758344020678">17.5</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="67.73657996679579">67.74</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_massiliensis_GCA_900459365"><span style="display:inline;">Streptococcus massiliensis GCA 900459365</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="68.67665699560021">68.68</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_acidominimus_GCA_001921825"><span style="display:inline;">Streptococcus acidominimus GCA 001921825</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="67.17651064879732">67.18</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_suis_NZ_LS483418"><span style="display:inline;">Streptococcus suis NZ LS483418</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="66.70646230836596">66.71</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_ovis_DSM_16829_GCA_000380125"><span style="display:inline;">Streptococcus ovis DSM 16829 GCA 000380125</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="66.34773071113148">66.35</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_orisratti_DSM_15617_GCA_000380105"><span style="display:inline;">Streptococcus orisratti DSM 15617 GCA 000380105</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="65.59783141149197">65.6</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_gallolyticus_GCA_900116655"><span style="display:inline;">Streptococcus gallolyticus GCA 900116655</span><sup title="assembly from synonym type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="65.52088907456479">65.52</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_minor_DSM_17118_GCA_000377005"><span style="display:inline;">Streptococcus minor DSM 17118 GCA 000377005</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="65.42807393828656">65.43</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_merionis_NZ_LT906439"><span style="display:inline;">Streptococcus merionis NZ LT906439</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="65.13568612164042">65.14</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_salivarius_NZ_LR134274"><span style="display:inline;">Streptococcus salivarius NZ LR134274</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="65.11872479693102">65.12</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_equinus_GCA_900459295"><span style="display:inline;">Streptococcus equinus GCA 900459295</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="65.09860909671062">65.1</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_infantarius_GCA_900459445"><span style="display:inline;">Streptococcus infantarius GCA 900459445</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.97269737312075">64.97</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_porci_DSM_23759_GCA_000423765"><span style="display:inline;">Streptococcus porci DSM 23759 GCA 000423765</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.87252534165506">64.87</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_ovuberis_GCA_012396585"><span style="display:inline;">Streptococcus ovuberis GCA 012396585</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.79772087874397">64.8</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_hyointestinalis_GCA_900459405"><span style="display:inline;">Streptococcus hyointestinalis GCA 900459405</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.73599496523782">64.74</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_pyogenes_NZ_CP040997"><span style="display:inline;">Streptococcus pyogenes NZ CP040997</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.53085153738517">64.53</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_agalactiae_GCA_900458965"><span style="display:inline;">Streptococcus agalactiae GCA 900458965</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.49787340136476">64.5</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_ratti_FA_1___DSM_20564_GCA_000286075"><span style="display:inline;">Streptococcus ratti FA 1   DSM 20564 GCA 000286075</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.36446372012448">64.36</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_pharyngis_GCA_007859195"><span style="display:inline;">Streptococcus pharyngis GCA 007859195</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.29185767869005">64.29</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_thoraltensis_DSM_12221_GCA_000380145"><span style="display:inline;">Streptococcus thoraltensis DSM 12221 GCA 000380145</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.21954852960565">64.22</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_ferus_DSM_20646_GCA_000372425"><span style="display:inline;">Streptococcus ferus DSM 20646 GCA 000372425</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.20114202505741">64.2</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_equi_GCA_900156215"><span style="display:inline;">Streptococcus equi GCA 900156215</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.11555084572078">64.12</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_urinalis_NZ_LR134323"><span style="display:inline;">Streptococcus urinalis NZ LR134323</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="64.06134405962987">64.06</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_criceti_HS_6_GCA_000187975"><span style="display:inline;">Streptococcus criceti HS 6 GCA 000187975</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.929042800451846">63.93</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_penaeicida_GCA_002887775"><span style="display:inline;">Streptococcus penaeicida GCA 002887775</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.760064278255456">63.76</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_sobrinus_DSM_20742___ATCC_33478_GCA_000439045"><span style="display:inline;">Streptococcus sobrinus DSM 20742   ATCC 33478 GCA 000439045</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.5927017109813">63.59</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_mutans_DSM_20523_GCA_000375505"><span style="display:inline;">Streptococcus mutans DSM 20523 GCA 000375505</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.55527986902792">63.56</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_pseudoporcinus_NZ_LR134341"><span style="display:inline;">Streptococcus pseudoporcinus NZ LR134341</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.552652555045064">63.55</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_halichoeri_GCA_009870755"><span style="display:inline;">Streptococcus halichoeri GCA 009870755</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.42083409971356">63.42</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_entericus_DSM_14446_GCA_000380025"><span style="display:inline;">Streptococcus entericus DSM 14446 GCA 000380025</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.41833681957384">63.42</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_macacae_NCTC_11558_GCA_900459485"><span style="display:inline;">Streptococcus macacae NCTC 11558 GCA 900459485</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.315102421065696">63.32</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_phocae_subsp__salmonis_GCA_000772915"><span style="display:inline;">Streptococcus phocae subsp. salmonis GCA 000772915</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.26336501527399">63.26</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#            <tr>
#                <td>
#                  <a href="http://microbial-genomes.org/projects/20/reference_datasets/Streptococcus_didelphis_DSM_15616_GCA_000380005"><span style="display:inline;">Streptococcus didelphis DSM 15616 GCA 000380005</span><sup title="assembly from type material" style="font-weight:bold;">T</sup></a>
#                </td>
#                <td>
#                  <span rel="tooltip" title="63.245642523433055">63.25</span>
#                </td>
#                <td>
#                  <span rel="tooltip" title="0.0">0.0</span>
#                </td>
#                <td>
#                  <span class="text-muted">(estimated)</span>
#                </td>
#              </tr>
#          </tbody>
#        </table>
#      </div>
#    </div>
#  </div>
#</div>
#
#</div>
#
#<div class="result-footer">
#  <i class="glyphicon glyphicon-download" title="downloads"> </i>
#  <a rel="tooltip" target="_blank" title="295 Bytes" href="http://microbial-genomes.org/query_datasets/42005/result/distances/aai_medoids">aai medoids</a>, <a rel="tooltip" target="_blank" title="19 KB" href="http://microbial-genomes.org/query_datasets/42005/result/distances/haai_db">haai db</a>, <a rel="tooltip" target="_blank" title="19 KB" href="http://microbial-genomes.org/query_datasets/42005/result/distances/aai_db">aai db</a>, <a rel="tooltip" target="_blank" title="2 KB" href="http://microbial-genomes.org/query_datasets/42005/result/distances/ani_db">ani db</a>, <a rel="tooltip" target="_blank" title="10.1 KB" href="http://microbial-genomes.org/query_datasets/42005/result/distances/ref_tree">ref tree</a>, <a rel="tooltip" target="_blank" title="45.5 KB" href="http://microbial-genomes.org/query_datasets/42005/result/distances/ref_tree_pdf">ref tree pdf</a>, or <a rel="tooltip" target="_blank" title="902 Bytes" href="http://microbial-genomes.org/query_datasets/42005/result/distances/intax_test">intax test</a>
#  <br>
#  <i class="glyphicon glyphicon-info-sign" title="downloads"> </i>
#  <a href="http://manual.microbial-genomes.org/part5/workflow.html#distances">Learn more</a>
#  <div class="timestamp">18 days ago</div>
#</div>
#</div>
#<div class="modal fade" id="info-msg-96387d7d-4377-4300-ac3b-b7014a337338" tabindex="-1" role="dialog" aria-labelledby="info-msg-96387d7d-4377-4300-ac3b-b7014a337338-h"><div class="modal-dialog " role="document"><div class="modal-content "><div class="modal-header"><button name="button" type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">×</span></button><h4 class="modal-title" id="info-msg-96387d7d-4377-4300-ac3b-b7014a337338-h">Taxonomic classification of query datasets</h4></div><div class="modal-body">
#    <p>
#      The taxonomic classification of your dataset is inferred by the maximum
#      Average Amino Acid Identity (AAI) found against all the genomes in the
#      database. The p-value is estimated from the empirical distribution
#      observed in all the reference genomes of NCBI's RefSeq at each taxonomic 
#      level, and indicates the probability of a different classification with
#      the observed AAI.
#    </p><p>
#      More formally, the p-value indicates the probability of
#      any two genomes not being classified in the same taxon at that taxonomic
#      rank, given that the AAI between them is greater than or equal to the one
#      observed. For example, a p-value of 0.01 for a genus classification
#      indicates that only 1% of all the pairs of genomes in RefSeq with AAI
#      above the observed value are from different genera; therefore the query
#      dataset most likely belongs to the same genus as the closest relative in
#      the database.
#    </p>
#  </div></div></div></div><div class="modal fade" id="info-msg-fa40afad-5f94-4847-aa94-f7ecbb85b004" tabindex="-1" role="dialog" aria-labelledby="info-msg-fa40afad-5f94-4847-aa94-f7ecbb85b004-h"><div class="modal-dialog " role="document"><div class="modal-content "><div class="modal-header"><button name="button" type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">×</span></button><h4 class="modal-title" id="info-msg-fa40afad-5f94-4847-aa94-f7ecbb85b004-h">Taxonomic novelty of query datasets</h4></div><div class="modal-body">
#    <p>
#      The taxonomic novelty of the query dataset is an analysis indicating the
#      taxonomic rank at which your dataset represents a novel taxon with respect
#      to the database. Note that the database does not contain genomes for all
#      described species, and the query dataset may belong to a previously
#      described taxon not represented in the database.
#    </p><p>
#      The taxonomic novelty is determined by the maximum Average Amino Acid
#      Identity (AAI) found against the genomes in the database. The p-value is
#      estimated from the empirical distribution observed in all the reference
#      genomes of NCBI's RefSeq at each taxonomic level, and indicates the
#      probability of the observed AAI between genomes in the same taxon.
#    </p><p>
#      More formally, the p-value indicates the probability of any two genomes
#      having an AAI lesser than or equal to the one observed, given that they
#      are classified in the same taxon at that taxonomic rank. For example, a
#      p-value of 0.01 for a novel family indicates that only 1% of all the
#      pairs of genomes in RefSeq classified in the same family have an AAI below
#      the observed value; therefore the query dataset most likely belongs to a
#      family not represented in the database.
#    </p>
#  </div></div></div></div><div class="modal fade" id="info-msg-0d90c62d-ffa8-40ce-9b18-a9609fbd64d4" tabindex="-1" role="dialog" aria-labelledby="info-msg-0d90c62d-ffa8-40ce-9b18-a9609fbd64d4-h"><div class="modal-dialog " role="document"><div class="modal-content "><div class="modal-header"><button name="button" type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">×</span></button><h4 class="modal-title" id="info-msg-0d90c62d-ffa8-40ce-9b18-a9609fbd64d4-h">Genome relatedness to database references</h4></div><div class="modal-body">
#  <p>
#    Average sequence identities to reference datasets in the database. Displays
#    only the top-50 values.
#  </p>
#  <p><b>ANI:</b>Average Nucleotide Identity</p>
#  <p><b>AAI:</b>Average Amino Acid Identity</p>
#</div></div></div></div>
#</div>
#        <div id="query_dataset_result_ssu_42005">
#
#<h3 role="button" data-toggle="collapse" class="result-title" id="result-ssu" href="#result-cnt-ssu" aria-expanded="false" aria-controls="result-cnt-ssu">
#  <i class="glyphicon glyphicon-tasks"> </i>
#  Ribosomal RNA (small subunit)
#</h3>
#
#<div class="result collapse" id="result-cnt-ssu" aria-expanded="false">
#<div class="result-body">
#        <p>MiGA didn't detect SSU sequences.</p>
#
#    <p>
#      <b>Ssu:</b>
#      0
#      <br>
#      <b>Complete ssu:</b>
#      0
#      <br>
#    </p>
#</div>
#
#<div class="result-footer">
#  <i class="glyphicon glyphicon-download" title="downloads"> </i>
#  <a rel="tooltip" target="_blank" title="84 Bytes" href="http://microbial-genomes.org/query_datasets/42005/result/ssu/gff">gff</a> or <a rel="tooltip" target="_blank" title="71 Bytes" href="http://microbial-genomes.org/query_datasets/42005/result/ssu/all_ssu_genes">all ssu genes</a>
#  <br>
#  <i class="glyphicon glyphicon-info-sign" title="downloads"> </i>
#  <a href="http://manual.microbial-genomes.org/part5/workflow.html#ssu">Learn more</a>
#  <div class="timestamp">about 1 month ago</div>
#</div>
#</div>
#
#</div>
#        <div id="query_dataset_result_essential_genes_42005">
#
#<h3 role="button" data-toggle="collapse" class="result-title" id="result-essential_genes" href="#result-cnt-essential_genes" aria-expanded="false" aria-controls="result-cnt-essential_genes">
#  <i class="glyphicon glyphicon-tasks"> </i>
#  Quality (essential genes)
#</h3>
#
#<div class="result collapse" id="result-cnt-essential_genes" aria-expanded="false">
#<div class="result-body">
#    <div style="margin: 2em 0;"><div class="row" style="margin-bottom:10px;"><strong class="col-sm-4 text-right">Completeness</strong><div class="col-sm-4"><div class="progress" style="margin-bottom:0;"><div class="progress-bar progress-bar-success" role="progressbar" aria-valuenow="100.0" aria-valuemin="0" aria-valuemax="100" style="width: 100.0%;"> </div></div></div><div class="col-sm-4 text-success">100.0% (very high)</div></div><div class="row" style="margin-bottom:10px;"><strong class="col-sm-4 text-right">Contamination</strong><div class="col-sm-4"><div class="progress" style="margin-bottom:0;"><div class="progress-bar progress-bar-success" role="progressbar" aria-valuenow="1.9" aria-valuemin="0" aria-valuemax="100" style="width: 1.9%;"> </div></div></div><div class="col-sm-4 text-success">1.9% (very low)</div></div><div class="row" style="margin-bottom:10px;"><hr><strong class="col-sm-4 text-right">Quality</strong><div class="col-sm-4"><div class="progress" style="margin-bottom:0;"><div class="progress-bar progress-bar-success" role="progressbar" aria-valuenow="90.5" aria-valuemin="0" aria-valuemax="100" style="width: 90.5%;"> </div></div></div><div class="col-sm-4 text-success">90.5% (excellent)</div></div></div>
#
#<p>Genes commonly found in Bacteria and Archaea detected:</p>
#  <div class="panel panel-default" id="essential_genes">
#    <div class="panel-heading" role="tab" id="essential_genes-h">
#      <h4 class="panel-title">
#        <a role="button" data-toggle="collapse" data-parent="#essential_genes" class="btn" href="#essential_genes-b" aria-expanded="false" aria-controls="essential_genes-b">
#          Full report
#        </a>
#      </h4>
#    </div>
#    <div id="essential_genes-b" class="panel-collapse collapse" role="tabpanel" aria-labelledby="essential_genes-h">
#      <div class="panel-body" style="padding-top:0;">
#        <br>
#             <b>Collection:</b> dupont_2012 B
	elsif($line=~m/<b>Collection:<\/b>/){
		$line=~s/.*<b>Collection:<\/b> //;
		print OUTqual "$line\t";
	}
#
#                <br>
#             <b>Essential genes found:</b> 106/106.
	elsif($line=~m/<b>Essential genes found:<\/b>/){
		$line=~s/.*<b>Essential genes found:<\/b> //;
		$line=~s/\.//;
		my @temp=split("/",$line);
		print OUTqual "$temp[0]\t$temp[1]\t";
	}

#
#                <br>
#             <b>Completeness:</b> 100.0%.
	elsif($line=~m/<b>Completeness:<\/b>/ && $comp==0){
		$line=~s/.*> //;
		$line=~s/%\.//;
		print OUTqual "$line\t";
		$comp=1;	#occurs twice
	}
#
#                <br>
#             <b>Contamination:</b> 1.9%.
	elsif($line=~m/<b>Contamination:<\/b>/ && $cont==0){
		$line=~s/.*> //;
		$line=~s/%\.//;
		print OUTqual "$line\t";
		$cont=1;	#occurs twice
	}
#
#                <br>
#             <b>Multiple copies:</b> 
#	elsif($line=~m/<b>Multiple copies:<\/b>/){
#		$line=~s/.*<b>Multiple copies:<\/b> //;
#		print OUTqual "$line\t";
#	}
#
#                <br>
#              &nbsp;&nbsp;<b>2 <a target="_blank" rel="noopener" href="http://pfam.xfam.org/family/Methyltransf_5">Methyltransf_5</a>: </b>
#                MraW methylase family.
#        <br>
#              &nbsp;&nbsp;<b>2 <a target="_blank" rel="noopener" href="http://www.jcvi.org/cgi-bin/tigrfams/HmmReportPage.cgi?acc=TIGR00436">TIGR00436</a>: era: </b>
#                GTP-binding protein Era.
#</div>
#    </div>
#  </div>
#
#    <p>
#      <b>Completeness:</b>
#      100 %
#      <br>
#      <b>Contamination:</b>
#      1.9 %
#      <br>
#      <b>Quality:</b>
	elsif($line=~m/<b>Quality:<\/b>/){
		$quality=$numlin+1;
	}
	elsif($quality==$numlin){
		$line=~s/ +//g;
		print OUTqual "$line\t";
	}

#      90.5
#      <br>
#    </p>
#</div>
#
#<div class="result-footer">
#  <i class="glyphicon glyphicon-download" title="downloads"> </i>
#  <a rel="tooltip" target="_blank" title="39.7 KB" href="http://microbial-genomes.org/query_datasets/42005/result/essential_genes/ess_genes">ess genes</a>, <a rel="tooltip" target="_blank" title="Folder" href="http://microbial-genomes.org/query_datasets/42005/result/essential_genes/collection">collection</a>, <a rel="tooltip" target="_blank" title="221 Bytes" href="http://microbial-genomes.org/query_datasets/42005/result/essential_genes/report">report</a>, <a rel="tooltip" target="_blank" title="66.6 KB" href="http://microbial-genomes.org/query_datasets/42005/result/essential_genes/alignments">alignments</a>, or <a rel="tooltip" target="_blank" title="529 Bytes" href="http://microbial-genomes.org/query_datasets/42005/result/essential_genes/raw_report">raw report</a>
#  <br>
#  <i class="glyphicon glyphicon-info-sign" title="downloads"> </i>
#  <a href="http://manual.microbial-genomes.org/part5/workflow.html#essential-genes">Learn more</a>
#  <div class="timestamp">about 1 month ago</div>
#</div>
#</div>
#
#</div>
#        <div id="query_dataset_result_cds_42005">
#
#<h3 role="button" data-toggle="collapse" class="result-title" id="result-cds" href="#result-cnt-cds" aria-expanded="false" aria-controls="result-cnt-cds">
#  <i class="glyphicon glyphicon-tasks"> </i>
#  Gene prediction
#</h3>
#
#<div class="result collapse" id="result-cnt-cds" aria-expanded="false">
#<div class="result-body">
#    <p>
#      <b>Predicted proteins:</b>
#      1,807
#      <br>
#      <b>Average length:</b>
#      313.2767 aa
#      <br>
#      <b>Coding density:</b>
#      89.537 %
#      <br>
#      <b>Codon table:</b>
#      11
#      <br>
#    </p>
#</div>
#
#<div class="result-footer">
#  <i class="glyphicon glyphicon-download" title="downloads"> </i>
#  <a rel="tooltip" target="_blank" title="335 KB" href="http://microbial-genomes.org/query_datasets/42005/result/cds/proteins">proteins</a>, <a rel="tooltip" target="_blank" title="499 KB" href="http://microbial-genomes.org/query_datasets/42005/result/cds/genes">genes</a>, or <a rel="tooltip" target="_blank" title="66.3 KB" href="http://microbial-genomes.org/query_datasets/42005/result/cds/gff3">gff3</a>
#  <br>
#  <i class="glyphicon glyphicon-info-sign" title="downloads"> </i>
#  <a href="http://manual.microbial-genomes.org/part5/workflow.html#cds">Learn more</a>
#  <div class="timestamp">about 1 month ago</div>
#</div>
#</div>
#
#</div>
#        <div id="query_dataset_result_assembly_42005">
#
#<h3 role="button" data-toggle="collapse" class="result-title" id="result-assembly" href="#result-cnt-assembly" aria-expanded="false" aria-controls="result-cnt-assembly">
#  <i class="glyphicon glyphicon-tasks"> </i>
#  Assembly
#</h3>
#
#<div class="result collapse" id="result-cnt-assembly" aria-expanded="false">
#<div class="result-body">
#    <div id="d58d8c3c-8d72-483f-acc5-d707bb962bde" class="js-plotly-plot"><div class="plot-container plotly"><div class="user-select-none svg-container" style="position: relative; width: 700px; height: 450px;"><svg class="main-svg" xmlns="http://www.w3.org/2000/svg" xlink="http://www.w3.org/1999/xlink" width="700" height="450" style="background: rgb(255, 255, 255) none repeat scroll 0% 0%;"><defs id="defs-6bf225"><g class="clips"><clipPath id="clip6bf225xyplot" class="plotclip"><rect width="540" height="270"></rect></clipPath><clipPath class="axesclip" id="clip6bf225x"><rect x="80" y="0" width="540" height="450"></rect></clipPath><clipPath class="axesclip" id="clip6bf225y"><rect x="0" y="100" width="700" height="270"></rect></clipPath><clipPath class="axesclip" id="clip6bf225xy"><rect x="80" y="100" width="540" height="270"></rect></clipPath></g><g class="gradients"></g></defs><g class="bglayer"></g><g class="draglayer cursor-crosshair"><g class="xy"><rect class="nsewdrag drag" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="80" y="100" width="540" height="270"></rect><rect class="nwdrag drag cursor-nw-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="60" y="80" width="20" height="20"></rect><rect class="nedrag drag cursor-ne-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="620" y="80" width="20" height="20"></rect><rect class="swdrag drag cursor-sw-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="60" y="370" width="20" height="20"></rect><rect class="sedrag drag cursor-se-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="620" y="370" width="20" height="20"></rect><rect class="ewdrag drag cursor-ew-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="134" y="370.5" width="432" height="20"></rect><rect class="wdrag drag cursor-w-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="80" y="370.5" width="54" height="20"></rect><rect class="edrag drag cursor-e-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="566" y="370.5" width="54" height="20"></rect><rect class="nsdrag drag cursor-ns-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="59.5" y="127" width="20" height="216"></rect><rect class="sdrag drag cursor-s-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="59.5" y="343" width="20" height="27"></rect><rect class="ndrag drag cursor-n-resize" style="fill: transparent; stroke-width: 0px; pointer-events: all;" data-subplot="xy" x="59.5" y="100" width="20" height="27"></rect></g></g><g class="layer-below"><g class="imagelayer"></g><g class="shapelayer"></g></g><g class="cartesianlayer"><g class="subplot xy"><g class="layer-subplot"><g class="shapelayer"></g><g class="imagelayer"></g></g><g class="gridlayer"><g class="x"><path class="xgrid crisp" transform="translate(95.61,0)" d="M0,100v270" style="stroke: rgb(238, 238, 238); stroke-opacity: 1; stroke-width: 1px;"></path><path class="xgrid crisp" transform="translate(171.93,0)" d="M0,100v270" style="stroke: rgb(238, 238, 238); stroke-opacity: 1; stroke-width: 1px;"></path><path class="xgrid crisp" transform="translate(248.26,0)" d="M0,100v270" style="stroke: rgb(238, 238, 238); stroke-opacity: 1; stroke-width: 1px;"></path><path class="xgrid crisp" transform="translate(324.58000000000004,0)" d="M0,100v270" style="stroke: rgb(238, 238, 238); stroke-opacity: 1; stroke-width: 1px;"></path><path class="xgrid crisp" transform="translate(400.9,0)" d="M0,100v270" style="stroke: rgb(238, 238, 238); stroke-opacity: 1; stroke-width: 1px;"></path><path class="xgrid crisp" transform="translate(477.23,0)" d="M0,100v270" style="stroke: rgb(238, 238, 238); stroke-opacity: 1; stroke-width: 1px;"></path><path class="xgrid crisp" transform="translate(553.55,0)" d="M0,100v270" style="stroke: rgb(238, 238, 238); stroke-opacity: 1; stroke-width: 1px;"></path></g><g class="y"></g></g><g class="zerolinelayer"></g><path class="xlines-below"></path><path class="ylines-below"></path><g class="overlines-below"></g><g class="xaxislayer-below"></g><g class="yaxislayer-below"></g><g class="overaxes-below"></g><g class="plot" transform="translate(80, 100)" clip-path="url('#clip6bf225xyplot')"><g class="boxlayer mlayer"><g class="trace boxes" style="opacity: 1;"><path style="vector-effect: non-scaling-stroke; stroke-width: 2px; stroke: rgb(31, 119, 180); stroke-opacity: 1; fill: rgb(31, 119, 180); fill-opacity: 0.5;" class="box" d="M145.44,201.15V68.85M82.59,201.15V68.85H273.16V201.15ZM82.59,135H27M273.16,135H513M27,168.08V101.93M513,168.08V101.93"></path><g class="points"></g></g></g></g><g class="overplot"></g><path class="xlines-above crisp" style="fill: none;" d="M0,0"></path><path class="ylines-above crisp" style="fill: none;" d="M0,0"></path><g class="overlines-above"></g><g class="xaxislayer-above"><g class="xtick"><text text-anchor="middle" x="0" y="383" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="0" data-math="N" transform="translate(95.61,0)">0</text></g><g class="xtick"><text text-anchor="middle" x="0" y="383" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="2M" data-math="N" transform="translate(171.93,0)">2M</text></g><g class="xtick"><text text-anchor="middle" x="0" y="383" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="4M" data-math="N" transform="translate(248.26,0)">4M</text></g><g class="xtick"><text text-anchor="middle" x="0" y="383" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="6M" data-math="N" transform="translate(324.58000000000004,0)">6M</text></g><g class="xtick"><text text-anchor="middle" x="0" y="383" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="8M" data-math="N" transform="translate(400.9,0)">8M</text></g><g class="xtick"><text text-anchor="middle" x="0" y="383" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="10M" data-math="N" transform="translate(477.23,0)">10M</text></g><g class="xtick"><text text-anchor="middle" x="0" y="383" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="12M" data-math="N" transform="translate(553.55,0)">12M</text></g></g><g class="yaxislayer-above"><g class="ytick"><text text-anchor="end" x="79" y="4.199999999999999" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" data-unformatted="RefSeq" data-math="N" transform="translate(0,235)">RefSeq</text></g></g><g class="overaxes-above"></g></g></g><g class="polarlayer"></g><g class="ternarylayer"></g><g class="geolayer"></g><g class="funnelarealayer"></g><g class="pielayer"></g><g class="treemaplayer"></g><g class="sunburstlayer"></g><g class="glimages"></g></svg><div class="gl-container"></div><svg class="main-svg" xmlns="http://www.w3.org/2000/svg" xlink="http://www.w3.org/1999/xlink" width="700" height="450"><defs id="topdefs-6bf225"><g class="clips"></g></defs><g class="indicatorlayer"></g><g class="layer-above"><g class="imagelayer"></g><g class="shapelayer"></g></g><g class="infolayer"><g class="g-gtitle"></g><g class="g-xtitle"><text class="xtitle" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 14px; fill: rgb(68, 68, 68); opacity: 1; font-weight: normal; white-space: pre;" x="350" y="422" text-anchor="middle" data-unformatted="Length (bp)" data-math="N">Length (bp)</text></g><g class="g-ytitle"></g><g class="annotation" data-index="0" style="opacity: 1;"><g class="annotation-text-g" transform="rotate(0,167.99,117.5)"><g class="cursor-pointer" transform="translate(134, 108)"><rect class="bg" style="stroke-width: 1px; stroke: rgb(0, 0, 0); stroke-opacity: 0; fill: rgb(0, 0, 0); fill-opacity: 0;" x="0.5" y="0.5" width="66" height="19"></rect><text class="annotation-text" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" text-anchor="middle" data-unformatted="Total length" data-math="N" x="33.6669921875" y="14">Total length</text></g></g><g style="opacity: 1;" class="annotation-arrow-g"><path d="M167.99,127.5L167.99,167.5" style="stroke-width: 2px; stroke: rgb(68, 68, 68); stroke-opacity: 1; stroke-dasharray: 0px, 0px, 38.8px, 40px;"></path><path d="M-2.4,-3V3L0.6,0Z" transform="translate(167.99000549316406,166.3000030517578)rotate(90)scale(2)" style="fill: rgb(68, 68, 68); stroke-width: 0px;"></path></g></g><g class="annotation" data-index="1" style="opacity: 1;"><g class="annotation-text-g" transform="rotate(0,100.61,137.5)"><g class="cursor-pointer" transform="translate(87, 128)"><rect class="bg" style="stroke-width: 1px; stroke: rgb(0, 0, 0); stroke-opacity: 0; fill: rgb(0, 0, 0); fill-opacity: 0;" x="0.5" y="0.5" width="27" height="19"></rect><text class="annotation-text" style="font-family: &quot;Open Sans&quot;, verdana, arial, sans-serif; font-size: 12px; fill: rgb(68, 68, 68); fill-opacity: 1; white-space: pre;" text-anchor="middle" data-unformatted="N50" data-math="N" x="14" y="14">N50</text></g></g><g style="opacity: 1;" class="annotation-arrow-g"><path d="M100.61,147.5L100.61,167.5" style="stroke-width: 2px; stroke: rgb(68, 68, 68); stroke-opacity: 1; stroke-dasharray: 0px, 0px, 18.8px, 20px;"></path><path d="M-2.4,-3V3L0.6,0Z" transform="translate(100.61000061035156,166.3000030517578)rotate(90)scale(2)" style="fill: rgb(68, 68, 68); stroke-width: 0px;"></path></g></g></g><g class="menulayer"></g><g class="zoomlayer"></g></svg><div class="modebar-container" style="position: absolute; top: 0px; right: 0px; width: 100%;"><div id="modebar-6bf225" class="modebar modebar--hover ease-bg"><div class="modebar-group"><a rel="tooltip" class="modebar-btn" data-title="Download plot as a png" data-toggle="false" data-gravity="n"><svg viewBox="0 0 1000 1000" class="icon" height="1em" width="1em"><path d="m500 450c-83 0-150-67-150-150 0-83 67-150 150-150 83 0 150 67 150 150 0 83-67 150-150 150z m400 150h-120c-16 0-34 13-39 29l-31 93c-6 15-23 28-40 28h-340c-16 0-34-13-39-28l-31-94c-6-15-23-28-40-28h-120c-55 0-100-45-100-100v-450c0-55 45-100 100-100h800c55 0 100 45 100 100v450c0 55-45 100-100 100z m-400-550c-138 0-250 112-250 250 0 138 112 250 250 250 138 0 250-112 250-250 0-138-112-250-250-250z m365 380c-19 0-35 16-35 35 0 19 16 35 35 35 19 0 35-16 35-35 0-19-16-35-35-35z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a></div><div class="modebar-group"><a rel="tooltip" class="modebar-btn active" data-title="Zoom" data-attr="dragmode" data-val="zoom" data-toggle="false" data-gravity="n"><svg viewBox="0 0 1000 1000" class="icon" height="1em" width="1em"><path d="m1000-25l-250 251c40 63 63 138 63 218 0 224-182 406-407 406-224 0-406-182-406-406s183-406 407-406c80 0 155 22 218 62l250-250 125 125z m-812 250l0 438 437 0 0-438-437 0z m62 375l313 0 0-312-313 0 0 312z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a><a rel="tooltip" class="modebar-btn" data-title="Pan" data-attr="dragmode" data-val="pan" data-toggle="false" data-gravity="n"><svg viewBox="0 0 1000 1000" class="icon" height="1em" width="1em"><path d="m1000 350l-187 188 0-125-250 0 0 250 125 0-188 187-187-187 125 0 0-250-250 0 0 125-188-188 186-187 0 125 252 0 0-250-125 0 187-188 188 188-125 0 0 250 250 0 0-126 187 188z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a></div><div class="modebar-group"><a rel="tooltip" class="modebar-btn" data-title="Zoom in" data-attr="zoom" data-val="in" data-toggle="false" data-gravity="n"><svg viewBox="0 0 875 1000" class="icon" height="1em" width="1em"><path d="m1 787l0-875 875 0 0 875-875 0z m687-500l-187 0 0-187-125 0 0 187-188 0 0 125 188 0 0 187 125 0 0-187 187 0 0-125z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a><a rel="tooltip" class="modebar-btn" data-title="Zoom out" data-attr="zoom" data-val="out" data-toggle="false" data-gravity="n"><svg viewBox="0 0 875 1000" class="icon" height="1em" width="1em"><path d="m0 788l0-876 875 0 0 876-875 0z m688-500l-500 0 0 125 500 0 0-125z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a><a rel="tooltip" class="modebar-btn" data-title="Autoscale" data-attr="zoom" data-val="auto" data-toggle="false" data-gravity="n"><svg viewBox="0 0 1000 1000" class="icon" height="1em" width="1em"><path d="m250 850l-187 0-63 0 0-62 0-188 63 0 0 188 187 0 0 62z m688 0l-188 0 0-62 188 0 0-188 62 0 0 188 0 62-62 0z m-875-938l0 188-63 0 0-188 0-62 63 0 187 0 0 62-187 0z m875 188l0-188-188 0 0-62 188 0 62 0 0 62 0 188-62 0z m-125 188l-1 0-93-94-156 156 156 156 92-93 2 0 0 250-250 0 0-2 93-92-156-156-156 156 94 92 0 2-250 0 0-250 0 0 93 93 157-156-157-156-93 94 0 0 0-250 250 0 0 0-94 93 156 157 156-157-93-93 0 0 250 0 0 250z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a><a rel="tooltip" class="modebar-btn" data-title="Reset axes" data-attr="zoom" data-val="reset" data-toggle="false" data-gravity="n"><svg viewBox="0 0 928.6 1000" class="icon" height="1em" width="1em"><path d="m786 296v-267q0-15-11-26t-25-10h-214v214h-143v-214h-214q-15 0-25 10t-11 26v267q0 1 0 2t0 2l321 264 321-264q1-1 1-4z m124 39l-34-41q-5-5-12-6h-2q-7 0-12 3l-386 322-386-322q-7-4-13-4-7 2-12 7l-35 41q-4 5-3 13t6 12l401 334q18 15 42 15t43-15l136-114v109q0 8 5 13t13 5h107q8 0 13-5t5-13v-227l122-102q5-5 6-12t-4-13z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a></div><div class="modebar-group"><a rel="tooltip" class="modebar-btn" data-title="Toggle Spike Lines" data-attr="_cartesianSpikesEnabled" data-val="on" data-toggle="false" data-gravity="n"><svg viewBox="0 0 1000 1000" class="icon" height="1em" width="1em"><path d="M512 409c0-57-46-104-103-104-57 0-104 47-104 104 0 57 47 103 104 103 57 0 103-46 103-103z m-327-39l92 0 0 92-92 0z m-185 0l92 0 0 92-92 0z m370-186l92 0 0 93-92 0z m0-184l92 0 0 92-92 0z" transform="matrix(1.5 0 0 -1.5 0 850)"></path></svg></a><a rel="tooltip" class="modebar-btn" data-title="Show closest data on hover" data-attr="hovermode" data-val="closest" data-toggle="false" data-gravity="ne"><svg viewBox="0 0 1500 1000" class="icon" height="1em" width="1em"><path d="m375 725l0 0-375-375 375-374 0-1 1125 0 0 750-1125 0z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a><a rel="tooltip" class="modebar-btn active" data-title="Compare data on hover" data-attr="hovermode" data-val="y" data-toggle="false" data-gravity="ne"><svg viewBox="0 0 1125 1000" class="icon" height="1em" width="1em"><path d="m187 786l0 2-187-188 188-187 0 0 937 0 0 373-938 0z m0-499l0 1-187-188 188-188 0 0 937 0 0 376-938-1z" transform="matrix(1 0 0 -1 0 850)"></path></svg></a></div><div class="modebar-group"><a href="https://plotly.com/" target="_blank" data-title="Produced with Plotly" class="modebar-btn plotlyjsicon modebar-btn--logo"><svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 132 132" height="1em" width="1em"><defs><style>.cls-1 {fill: #3f4f75;} .cls-2 {fill: #80cfbe;} .cls-3 {fill: #fff;}</style></defs><title>plotly-logomark</title><g id="symbol"><rect class="cls-1" width="132" height="132" rx="6" ry="6"></rect><circle class="cls-2" cx="78" cy="54" r="6"></circle><circle class="cls-2" cx="102" cy="30" r="6"></circle><circle class="cls-2" cx="78" cy="30" r="6"></circle><circle class="cls-2" cx="54" cy="30" r="6"></circle><circle class="cls-2" cx="30" cy="30" r="6"></circle><circle class="cls-2" cx="30" cy="54" r="6"></circle><path class="cls-3" d="M30,72a6,6,0,0,0-6,6v24a6,6,0,0,0,12,0V78A6,6,0,0,0,30,72Z"></path><path class="cls-3" d="M78,72a6,6,0,0,0-6,6v24a6,6,0,0,0,12,0V78A6,6,0,0,0,78,72Z"></path><path class="cls-3" d="M54,48a6,6,0,0,0-6,6v48a6,6,0,0,0,12,0V54A6,6,0,0,0,54,48Z"></path><path class="cls-3" d="M102,48a6,6,0,0,0-6,6v48a6,6,0,0,0,12,0V54A6,6,0,0,0,102,48Z"></path></g></svg></a></div></div></div><svg class="main-svg" xmlns="http://www.w3.org/2000/svg" xlink="http://www.w3.org/1999/xlink" width="700" height="450"><g class="hoverlayer"></g></svg></div></div><script>
#//<![CDATA[
#$( Plotly.plot('d58d8c3c-8d72-483f-acc5-d707bb962bde', [{"x":[298471,2240716,3402093,4653910,13033779],"type":"box","hoverinfo":"none","name":"RefSeq"}], {"showlegend":false,"xaxis":{"autorange":true,"zeroline":false,"title":"Length (bp)"},"annotations":[{"x":1896728,"y":0.25,"text":"Total length","showarrow":true,"yanchor":"bottom","ax":0,"ay":-40},{"x":131071,"y":0.25,"text":"N50","showarrow":true,"yanchor":"bottom","ax":0,"ay":-20}]}) );
#//]]>
#</script></div>
#    <p>
#      <b>Contigs:</b>
#      33
#      <br>
#      <b>N50:</b>
#      131,071 bp
#      <br>
#      <b>Total length:</b>
#      1,896,728 bp
#      <br>
#      <b>Longest sequence:</b>
#      381,131 bp
#      <br>
#      <b>G+C content:</b>
#      40.3249 %
#      <br>
#      <b>X content:</b>
#      0 %
#      <br>
#      <b>G-C skew:</b>
#      -3.2041 %
#      <br>
#      <b>A-T skew:</b>
#      -1.0636 %
#      <br>
#    </p>
#</div>
#
#<div class="result-footer">
#  <i class="glyphicon glyphicon-download" title="downloads"> </i>
#  <a rel="tooltip" target="_blank" title="1.84 MB" href="http://microbial-genomes.org/query_datasets/42005/result/assembly/largecontigs">largecontigs</a>
#  <br>
#  <i class="glyphicon glyphicon-info-sign" title="downloads"> </i>
#  <a href="http://manual.microbial-genomes.org/part5/workflow.html#assembly">Learn more</a>
#  <div class="timestamp">about 1 month ago</div>
#</div>
#</div>
#
#</div>
#    <script>
#        $.ajax(
#          { url: "/projects/20/result/stats?q_ds=42005" }
#        ).done(function(data){
#          $("#query_dataset_result_stats_42005").
#            html(data) });
#        $.ajax(
#          { url: "/projects/20/result/distances?q_ds=42005" }
#        ).done(function(data){
#          $("#query_dataset_result_distances_42005").
#            html(data) });
#        $.ajax(
#          { url: "/projects/20/result/ssu?q_ds=42005" }
#        ).done(function(data){
#          $("#query_dataset_result_ssu_42005").
#            html(data) });
#        $.ajax(
#          { url: "/projects/20/result/essential_genes?q_ds=42005" }
#        ).done(function(data){
#          $("#query_dataset_result_essential_genes_42005").
#            html(data) });
#        $.ajax(
#          { url: "/projects/20/result/cds?q_ds=42005" }
#        ).done(function(data){
#          $("#query_dataset_result_cds_42005").
#            html(data) });
#        $.ajax(
#          { url: "/projects/20/result/assembly?q_ds=42005" }
#        ).done(function(data){
#          $("#query_dataset_result_assembly_42005").
#            html(data) });
#    </script>
#      <br><br>
#      <a type="button" class="btn btn-default btn-block" href="http://microbial-genomes.org/query_datasets/42005/run_mytaxa_scan">Execute MyTaxa scan analysis</a>
#  </div>
#</div>
#
#    <div class="push"></div>
#  </div>
#  
#  <footer class="footer container-fluid">
#  <div class="row">
#    <div class="col-md-8">
#      <nav>
#        <ul class="left-col">
#          <li><a target="_blank" rel="noopener" href="https://help.microbial-genomes.org/">Help</a></li>
#          <li><a target="_blank" rel="noopener" href="http://miga.microbial-genomes.org/">News</a></li>
#          <li><a target="_blank" rel="noopener" href="http://support.microbial-genomes.org/">Support</a></li>
#          <li><a href="http://microbial-genomes.org/about">About</a></li>
#          <li><a target="_blank" rel="noopener" href="http://roadmap.microbial-genomes.org/">Roadmap</a></li>
#        </ul>
#      </nav>
#    </div>
#    <div class="col-md-4 text-right right-col">
#      <small>
#        <p>
#          <a target="_blank" rel="noopener" href="http://code.microbial-genomes.org/miga-web">MiGA Web by Luis M. Rodriguez-R</a><br>
#          Powered by <a target="_blank" rel="noopener" href="http://code.microbial-genomes.org/miga">MiGA 0.7.15.2 - lithograph</a>
#        </p>
#      </small>
#    </div>
#  </div>
#  <div class="row center">
#    <div class="col-sm-4">
#      <a target="_blank" rel="noopener" href="http://miga.microbial-genomes.org/"><img alt="MiGA logo" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/MiGA-01f27acf6cf8127cb73f92672b96dd6b31ab571eff01c282cead78d.png"></a>
#    </div>
#    <div class="col-sm-4">
#      <h4>Funded by:</h4>
#      <a target="_blank" rel="noopener" href="http://nsf.gov/">US National Science Foundation</a><br>
#      <a target="_blank" rel="noopener" href="http://nsf.gov/awardsearch/showAward?AWD_ID=1356288">Award #1356288</a>
#    </div>
#    <div class="col-sm-4">
#      <h4>A collaboration between:</h4>
#      <a target="_blank" rel="noopener" href="http://enve-omics.gatech.edu/"><img alt="Georgia Tech logo" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/GT-766fd6d446be3f64d1c59841ae72781b8241cfc7656ff59e65251d07b.png"></a>
#      <a target="_blank" rel="noopener" href="http://rdp.cme.msu.edu/"><img alt="RDP logo" src="Streptococcus%20mitis%2027098%208%2093%20|%20MiGA%20Online_files/RDP-83c5f63a6637b13e0798643603802638076bfa118c536898d1492529.png"></a>
#    </div>
#  </div>
#</footer>
#
#  
#  <script>
#    $(function(){ $("[rel=tooltip]").tooltip({ placement: 'top'}); });
#  </script>
#
#
#<svg id="js-plotly-tester" xmlns="http://www.w3.org/2000/svg" xlink="http://www.w3.org/1999/xlink" style="position: absolute; left: -10000px; top: -10000px; width: 9000px; height: 9000px; z-index: 1;"><path class="js-reference-point" d="M0,0H1V1H0Z" style="stroke-width: 0px; fill: black;"></path></svg></body></html>
	}
	print OUTqual "\n";
}
close OUTtax ;
close OUTqual ;
close OUTnovel ;
close OUTtaxpval ;
close OUTanovel;
