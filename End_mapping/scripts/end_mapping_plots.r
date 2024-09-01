library (ggplot2)
library(ggh4x)

plastid_starts = read.table("plots/plastid_starts.txt", header=TRUE)

mitochondrial_starts = read.table("plots/mitochondrial_starts.txt", header=TRUE)

nuclear_starts = read.table("plots/nuclear_starts.txt", header=TRUE)

plastid_ends = read.table("plots/plastid_ends.txt", header=TRUE)

mitochondrial_ends = read.table("plots/mitochondrial_ends.txt", header=TRUE)

nuclear_ends = read.table("plots/nuclear_ends.txt", header=TRUE)

plastid_start_plot = ggplot(data=plastid_starts, aes(x=SprinzlPos, y=TPM/1000, fill=Type)) + geom_col() + facet_wrap(~IsoDecoder, scales="free_y", ncol = 6) + theme_bw() + ylab("Read Count") +xlab ("Read Mapping Start Position") + scale_fill_manual(values=c("gray20","firebrick")) + theme(legend.position = "top", strip.text = element_text(size=5, margin = margin(2,2,2,2)), axis.text = element_text(size=5), legend.title = element_blank(), legend.text=element_text(size=7), axis.title = element_text(size=8,face="bold")) + facetted_pos_scales( y = list (IsoDecoder == "AlaTGC" ~ scale_y_continuous(limits = c(-94.6332905915606,94.6332905915606)),IsoDecoder == "ArgACG" ~ scale_y_continuous(limits = c(-43.6966040069682,43.6966040069682)),IsoDecoder == "ArgTCT" ~ scale_y_continuous(limits = c(-7.3828136182637,7.3828136182637)),IsoDecoder == "AsnGTT" ~ scale_y_continuous(limits = c(-57.9806853752168,57.9806853752168)),IsoDecoder == "AspGTC" ~ scale_y_continuous(limits = c(-37.6108793390696,37.6108793390696)),IsoDecoder == "CysGCA" ~ scale_y_continuous(limits = c(-19.914362037146,19.914362037146)),IsoDecoder == "GlnTTG" ~ scale_y_continuous(limits = c(-25.0717989441303,25.0717989441303)),IsoDecoder == "GluTTC" ~ scale_y_continuous(limits = c(-125.976925108859,125.976925108859)),IsoDecoder == "GlyGCC" ~ scale_y_continuous(limits = c(-10.9642885612853,10.9642885612853)),IsoDecoder == "GlyTCC" ~ scale_y_continuous(limits = c(-9.48224110378493,9.48224110378493)),IsoDecoder == "HisGTG" ~ scale_y_continuous(limits = c(-55.6520762148008,55.6520762148008)),IsoDecoder == "IleCAT" ~ scale_y_continuous(limits = c(-18.913707846648,18.913707846648)),IsoDecoder == "IleGAT" ~ scale_y_continuous(limits = c(-56.0270153591925,56.0270153591925)),IsoDecoder == "LeuCAA" ~ scale_y_continuous(limits = c(-26.2201924138013,26.2201924138013)),IsoDecoder == "LeuTAA" ~ scale_y_continuous(limits = c(-11.7849770335339,11.7849770335339)),IsoDecoder == "LeuTAG" ~ scale_y_continuous(limits = c(-30.4393383137458,30.4393383137458)),IsoDecoder == "LysTTT" ~ scale_y_continuous(limits = c(-6.45550919200238,6.45550919200238)),IsoDecoder == "Met(e)CAT" ~ scale_y_continuous(limits = c(-11.5899660653072,11.5899660653072)),IsoDecoder == "Met(i)CAT" ~ scale_y_continuous(limits = c(-11.480759489395,11.480759489395)),IsoDecoder == "PheGAA" ~ scale_y_continuous(limits = c(-62.5327481659006,62.5327481659006)),IsoDecoder == "ProTGG" ~ scale_y_continuous(limits = c(-5.76674772026781,5.76674772026781)),IsoDecoder == "SerGCT" ~ scale_y_continuous(limits = c(-5.2092753867597,5.2092753867597)),IsoDecoder == "SerGGA" ~ scale_y_continuous(limits = c(-16.9713102630297,16.9713102630297)),IsoDecoder == "SerTGA" ~ scale_y_continuous(limits = c(-17.0627528264711,17.0627528264711)),IsoDecoder == "ThrGGT" ~ scale_y_continuous(limits = c(-13.9935213245771,13.9935213245771)),IsoDecoder == "ThrTGT" ~ scale_y_continuous(limits = c(-10.1031820486465,10.1031820486465)),IsoDecoder == "TrpCCA" ~ scale_y_continuous(limits = c(-15.8964961010438,15.8964961010438)),IsoDecoder == "TyrGTA" ~ scale_y_continuous(limits = c(-36.207128667473,36.207128667473)),IsoDecoder == "ValGAC" ~ scale_y_continuous(limits = c(-2.3931353503193,2.3931353503193)),IsoDecoder == "ValTAC" ~ scale_y_continuous(limits = c(-14.0634050467744,14.0634050467744))))


mitochondrial_start_plot = ggplot(data=mitochondrial_starts, aes(x=SprinzlPos, y=TPM/1000, fill=Type)) + geom_col() + facet_wrap(~IsoDecoder, scales="free_y", ncol = 6) + theme_bw() + ylab("Read Count") +xlab ("Read Mapping Start Position") + scale_fill_manual(values=c("gray20","firebrick")) + theme(legend.position = "top", strip.text = element_text(size=5, margin = margin(2,2,2,2)), axis.text = element_text(size=5), legend.title = element_blank(), legend.text=element_text(size=7), axis.title = element_text(size=8,face="bold")) + facetted_pos_scales( y = list (IsoDecoder == "AsnGTT" ~ scale_y_continuous(limits = c(-147.626985445917,147.626985445917)),IsoDecoder == "CysGCA" ~ scale_y_continuous(limits = c(-448.127581620769,448.127581620769)),IsoDecoder == "GlnTTG" ~ scale_y_continuous(limits = c(-5.54787998031688,5.54787998031688)),IsoDecoder == "GluTTC" ~ scale_y_continuous(limits = c(-72.1273347368769,72.1273347368769)),IsoDecoder == "GlyGCC" ~ scale_y_continuous(limits = c(-70.3816097145652,70.3816097145652)),IsoDecoder == "HisGTG" ~ scale_y_continuous(limits = c(-89.0360460873362,89.0360460873362)),IsoDecoder == "IleCAT" ~ scale_y_continuous(limits = c(-11.9452273662618,11.9452273662618)),IsoDecoder == "LysTTT" ~ scale_y_continuous(limits = c(-28.2525443449033,28.2525443449033)),IsoDecoder == "Met(e)CAT" ~ scale_y_continuous(limits = c(-8.09560805111614,8.09560805111614)),IsoDecoder == "Met(i)CAT" ~ scale_y_continuous(limits = c(-4.0439226807083,4.0439226807083)),IsoDecoder == "ProTGG" ~ scale_y_continuous(limits = c(-35.4162308606246,35.4162308606246)),IsoDecoder == "SerGCT" ~ scale_y_continuous(limits = c(-25.547683187262,25.547683187262)),IsoDecoder == "SerTGA" ~ scale_y_continuous(limits = c(-21.4716393477239,21.4716393477239)),IsoDecoder == "TyrGTA" ~ scale_y_continuous(limits = c(-76.1230269704805,76.1230269704805))))


nuclear_start_plot = ggplot(data=nuclear_starts, aes(x=SprinzlPos, y=TPM/1000, fill=Type)) + geom_col() + facet_wrap(~IsoDecoder, scales="free_y", ncol = 6) + theme_bw() + ylab("Read Count") +xlab ("Read Mapping Start Position") + scale_fill_manual(values=c("gray20","firebrick")) + theme(legend.position = "top", strip.text = element_text(size=5, margin = margin(2,2,2,2)), axis.text = element_text(size=5), legend.title = element_blank(), legend.text=element_text(size=7), axis.title = element_text(size=8,face="bold")) + facetted_pos_scales( y = list (IsoDecoder == "AlaAGC" ~ scale_y_continuous(limits = c(-49.4911,49.4911)),IsoDecoder == "AlaCGC" ~ scale_y_continuous(limits = c(-17.6163,17.6163)),IsoDecoder == "AlaTGC" ~ scale_y_continuous(limits = c(-34.933,34.933)),IsoDecoder == "ArgACG" ~ scale_y_continuous(limits = c(-11.4342,11.4342)),IsoDecoder == "ArgCCG" ~ scale_y_continuous(limits = c(-6.5698,6.5698)),IsoDecoder == "ArgCCT" ~ scale_y_continuous(limits = c(-12.5034,12.5034)),IsoDecoder == "ArgTCG" ~ scale_y_continuous(limits = c(-14.9079,14.9079)),IsoDecoder == "ArgTCT" ~ scale_y_continuous(limits = c(-9.1987,9.1987)),IsoDecoder == "AsnGTT" ~ scale_y_continuous(limits = c(-15.513,15.513)),IsoDecoder == "AspGTC" ~ scale_y_continuous(limits = c(-44.2655,44.2655)),IsoDecoder == "CysGCA" ~ scale_y_continuous(limits = c(-1.6003,1.6003)),IsoDecoder == "GlnCTG" ~ scale_y_continuous(limits = c(-7.462,7.462)),IsoDecoder == "GlnTTG" ~ scale_y_continuous(limits = c(-6.3728,6.3728)),IsoDecoder == "GluCTC" ~ scale_y_continuous(limits = c(-5.8672,5.8672)),IsoDecoder == "GluTTC" ~ scale_y_continuous(limits = c(-11.9906,11.9906)),IsoDecoder == "GlyCCC" ~ scale_y_continuous(limits = c(-12.5642,12.5642)),IsoDecoder == "GlyGCC" ~ scale_y_continuous(limits = c(-27.6392,27.6392)),IsoDecoder == "GlyTCC" ~ scale_y_continuous(limits = c(-13.196,13.196)),IsoDecoder == "HisGTG" ~ scale_y_continuous(limits = c(-5.388,5.388)),IsoDecoder == "IleAAT" ~ scale_y_continuous(limits = c(-39.8689,39.8689)),IsoDecoder == "IleTAT" ~ scale_y_continuous(limits = c(-0.8931,0.8931)),IsoDecoder == "LeuAAG" ~ scale_y_continuous(limits = c(-31.8511,31.8511)),IsoDecoder == "LeuCAA" ~ scale_y_continuous(limits = c(-7.3532,7.3532)),IsoDecoder == "LeuCAG" ~ scale_y_continuous(limits = c(-6.1163,6.1163)),IsoDecoder == "LeuTAA" ~ scale_y_continuous(limits = c(-4.25,4.25)),IsoDecoder == "LeuTAG" ~ scale_y_continuous(limits = c(-28.8755,28.8755)),IsoDecoder == "LysCTT" ~ scale_y_continuous(limits = c(-18.2483,18.2483)),IsoDecoder == "LysTTT" ~ scale_y_continuous(limits = c(-22.8355,22.8355)),IsoDecoder == "Met(e)CAT" ~ scale_y_continuous(limits = c(-9.3517,9.3517)),IsoDecoder == "Met(i)CAT" ~ scale_y_continuous(limits = c(-0.2831,0.2831)),IsoDecoder == "PheGAA" ~ scale_y_continuous(limits = c(-15.3957,15.3957)),IsoDecoder == "ProAGG" ~ scale_y_continuous(limits = c(-4.4986,4.4986)),IsoDecoder == "ProCGG" ~ scale_y_continuous(limits = c(-4.3491,4.3491)),IsoDecoder == "ProTGG" ~ scale_y_continuous(limits = c(-6.0825,6.0825)),IsoDecoder == "SerAGA" ~ scale_y_continuous(limits = c(-10.7824,10.7824)),IsoDecoder == "SerCGA" ~ scale_y_continuous(limits = c(-3.4373,3.4373)),IsoDecoder == "SerGCT" ~ scale_y_continuous(limits = c(-5.0691,5.0691)),IsoDecoder == "SerTGA" ~ scale_y_continuous(limits = c(-3.3347,3.3347)),IsoDecoder == "ThrAGT" ~ scale_y_continuous(limits = c(-5.5824,5.5824)),IsoDecoder == "ThrCGT" ~ scale_y_continuous(limits = c(-7.8202,7.8202)),IsoDecoder == "ThrTGT" ~ scale_y_continuous(limits = c(-12.8464,12.8464)),IsoDecoder == "TrpCCA" ~ scale_y_continuous(limits = c(-175.312,175.312)),IsoDecoder == "TyrGTA" ~ scale_y_continuous(limits = c(-22.6675,22.6675)),IsoDecoder == "ValAAC" ~ scale_y_continuous(limits = c(-9.4499,9.4499)),IsoDecoder == "ValCAC" ~ scale_y_continuous(limits = c(-2.7288,2.7288)),IsoDecoder == "ValTAC" ~ scale_y_continuous(limits = c(-10.2779,10.2779))))


plastid_end_plot = ggplot(data=plastid_ends, aes(x=SprinzlPos, y=TPM/1000, fill=Type)) + geom_col() + facet_wrap(~IsoDecoder, scales="free_y", ncol = 6) + theme_bw() + ylab("Read Count") +xlab ("Read Mapping End Position") + scale_fill_manual(values=c("gray20","firebrick")) + theme(legend.position = "top", strip.text = element_text(size=5, margin = margin(2,2,2,2)), axis.text = element_text(size=5), legend.title = element_blank(), legend.text=element_text(size=7), axis.title = element_text(size=8,face="bold")) + facetted_pos_scales( y = list (IsoDecoder == "AlaTGC" ~ scale_y_continuous(limits = c(-65.881966981085,65.881966981085)),IsoDecoder == "ArgACG" ~ scale_y_continuous(limits = c(-46.4980754717199,46.4980754717199)),IsoDecoder == "ArgTCT" ~ scale_y_continuous(limits = c(-10.1010645371904,10.1010645371904)),IsoDecoder == "AsnGTT" ~ scale_y_continuous(limits = c(-82.9925798937233,82.9925798937233)),IsoDecoder == "AspGTC" ~ scale_y_continuous(limits = c(-20.2370561430661,20.2370561430661)),IsoDecoder == "CysGCA" ~ scale_y_continuous(limits = c(-21.5975378833337,21.5975378833337)),IsoDecoder == "GlnTTG" ~ scale_y_continuous(limits = c(-16.9583301293017,16.9583301293017)),IsoDecoder == "GluTTC" ~ scale_y_continuous(limits = c(-98.4618418641124,98.4618418641124)),IsoDecoder == "GlyGCC" ~ scale_y_continuous(limits = c(-7.83444262714666,7.83444262714666)),IsoDecoder == "GlyTCC" ~ scale_y_continuous(limits = c(-8.19000150745394,8.19000150745394)),IsoDecoder == "HisGTG" ~ scale_y_continuous(limits = c(-38.3255260596303,38.3255260596303)),IsoDecoder == "IleCAT" ~ scale_y_continuous(limits = c(-28.2943883229571,28.2943883229571)),IsoDecoder == "IleGAT" ~ scale_y_continuous(limits = c(-47.2880968984359,47.2880968984359)),IsoDecoder == "LeuCAA" ~ scale_y_continuous(limits = c(-20.8220245113959,20.8220245113959)),IsoDecoder == "LeuTAA" ~ scale_y_continuous(limits = c(-11.2653706948883,11.2653706948883)),IsoDecoder == "LeuTAG" ~ scale_y_continuous(limits = c(-26.033198943372,26.033198943372)),IsoDecoder == "LysTTT" ~ scale_y_continuous(limits = c(-7.58001676306179,7.58001676306179)),IsoDecoder == "Met(e)CAT" ~ scale_y_continuous(limits = c(-10.1786816320923,10.1786816320923)),IsoDecoder == "Met(i)CAT" ~ scale_y_continuous(limits = c(-6.7073169180506,6.7073169180506)),IsoDecoder == "PheGAA" ~ scale_y_continuous(limits = c(-49.7799519088127,49.7799519088127)),IsoDecoder == "ProTGG" ~ scale_y_continuous(limits = c(-4.88395701278198,4.88395701278198)),IsoDecoder == "SerGCT" ~ scale_y_continuous(limits = c(-3.75617753373766,3.75617753373766)),IsoDecoder == "SerGGA" ~ scale_y_continuous(limits = c(-9.26944181423515,9.26944181423515)),IsoDecoder == "SerTGA" ~ scale_y_continuous(limits = c(-14.8294079801755,14.8294079801755)),IsoDecoder == "ThrGGT" ~ scale_y_continuous(limits = c(-11.7687253871072,11.7687253871072)),IsoDecoder == "ThrTGT" ~ scale_y_continuous(limits = c(-8.59464159373864,8.59464159373864)),IsoDecoder == "TrpCCA" ~ scale_y_continuous(limits = c(-15.6904633575577,15.6904633575577)),IsoDecoder == "TyrGTA" ~ scale_y_continuous(limits = c(-41.8261598499756,41.8261598499756)),IsoDecoder == "ValGAC" ~ scale_y_continuous(limits = c(-2.43826562510361,2.43826562510361)),IsoDecoder == "ValTAC" ~ scale_y_continuous(limits = c(-10.5886738212777,10.5886738212777))))


mitochondrial_end_plot = ggplot(data=mitochondrial_ends, aes(x=SprinzlPos, y=TPM/1000, fill=Type)) + geom_col() + facet_wrap(~IsoDecoder, scales="free_y", ncol = 6) + theme_bw() + ylab("Read Count") +xlab ("Read Mapping End Position") + scale_fill_manual(values=c("gray20","firebrick")) + theme(legend.position = "top", strip.text = element_text(size=5, margin = margin(2,2,2,2)), axis.text = element_text(size=5), legend.title = element_blank(), legend.text=element_text(size=7), axis.title = element_text(size=8,face="bold")) + facetted_pos_scales( y = list (IsoDecoder == "AsnGTT" ~ scale_y_continuous(limits = c(-130.633713749411,130.633713749411)),IsoDecoder == "CysGCA" ~ scale_y_continuous(limits = c(-358.473949638953,358.473949638953)),IsoDecoder == "GlnTTG" ~ scale_y_continuous(limits = c(-4.63809185899445,4.63809185899445)),IsoDecoder == "GluTTC" ~ scale_y_continuous(limits = c(-50.1343193049821,50.1343193049821)),IsoDecoder == "GlyGCC" ~ scale_y_continuous(limits = c(-56.416720864359,56.416720864359)),IsoDecoder == "HisGTG" ~ scale_y_continuous(limits = c(-40.4244519441228,40.4244519441228)),IsoDecoder == "IleCAT" ~ scale_y_continuous(limits = c(-9.97230035471534,9.97230035471534)),IsoDecoder == "LysTTT" ~ scale_y_continuous(limits = c(-23.1313381385626,23.1313381385626)),IsoDecoder == "Met(e)CAT" ~ scale_y_continuous(limits = c(-8.41493191903198,8.41493191903198)),IsoDecoder == "Met(i)CAT" ~ scale_y_continuous(limits = c(-6.78191602442365,6.78191602442365)),IsoDecoder == "ProTGG" ~ scale_y_continuous(limits = c(-26.717772554079,26.717772554079)),IsoDecoder == "SerGCT" ~ scale_y_continuous(limits = c(-12.7862157210675,12.7862157210675)),IsoDecoder == "SerTGA" ~ scale_y_continuous(limits = c(-15.1096855451692,15.1096855451692)),IsoDecoder == "TyrGTA" ~ scale_y_continuous(limits = c(-48.7614323702965,48.7614323702965))))


nuclear_end_plot = ggplot(data=nuclear_ends, aes(x=SprinzlPos, y=TPM/1000, fill=Type)) + geom_col() + facet_wrap(~IsoDecoder, scales="free_y", ncol = 6) + theme_bw() + ylab("Read Count") +xlab ("Read Mapping End Position") + scale_fill_manual(values=c("gray20","firebrick")) + theme(legend.position = "top", strip.text = element_text(size=5, margin = margin(2,2,2,2)), axis.text = element_text(size=5), legend.title = element_blank(), legend.text=element_text(size=7), axis.title = element_text(size=8,face="bold")) + facetted_pos_scales( y = list (IsoDecoder == "AlaAGC" ~ scale_y_continuous(limits = c(-47.1449,47.1449)),IsoDecoder == "AlaCGC" ~ scale_y_continuous(limits = c(-17.0725,17.0725)),IsoDecoder == "AlaTGC" ~ scale_y_continuous(limits = c(-34.5381,34.5381)),IsoDecoder == "ArgACG" ~ scale_y_continuous(limits = c(-12.3266,12.3266)),IsoDecoder == "ArgCCG" ~ scale_y_continuous(limits = c(-8.2021,8.2021)),IsoDecoder == "ArgCCT" ~ scale_y_continuous(limits = c(-8.6611,8.6611)),IsoDecoder == "ArgTCG" ~ scale_y_continuous(limits = c(-14.8173,14.8173)),IsoDecoder == "ArgTCT" ~ scale_y_continuous(limits = c(-9.5589,9.5589)),IsoDecoder == "AsnGTT" ~ scale_y_continuous(limits = c(-11.1405,11.1405)),IsoDecoder == "AspGTC" ~ scale_y_continuous(limits = c(-15.827,15.827)),IsoDecoder == "CysGCA" ~ scale_y_continuous(limits = c(-1.8678,1.8678)),IsoDecoder == "GlnCTG" ~ scale_y_continuous(limits = c(-7.0563,7.0563)),IsoDecoder == "GlnTTG" ~ scale_y_continuous(limits = c(-3.4003,3.4003)),IsoDecoder == "GluCTC" ~ scale_y_continuous(limits = c(-3.3419,3.3419)),IsoDecoder == "GluTTC" ~ scale_y_continuous(limits = c(-9.1599,9.1599)),IsoDecoder == "GlyCCC" ~ scale_y_continuous(limits = c(-7.1107,7.1107)),IsoDecoder == "GlyGCC" ~ scale_y_continuous(limits = c(-24.2223,24.2223)),IsoDecoder == "GlyTCC" ~ scale_y_continuous(limits = c(-8.4889,8.4889)),IsoDecoder == "HisGTG" ~ scale_y_continuous(limits = c(-3.5523,3.5523)),IsoDecoder == "IleAAT" ~ scale_y_continuous(limits = c(-37.6673,37.6673)),IsoDecoder == "IleTAT" ~ scale_y_continuous(limits = c(-1.6361,1.6361)),IsoDecoder == "LeuAAG" ~ scale_y_continuous(limits = c(-17.5692,17.5692)),IsoDecoder == "LeuCAA" ~ scale_y_continuous(limits = c(-8.3803,8.3803)),IsoDecoder == "LeuCAG" ~ scale_y_continuous(limits = c(-3.3211,3.3211)),IsoDecoder == "LeuTAA" ~ scale_y_continuous(limits = c(-4.3542,4.3542)),IsoDecoder == "LeuTAG" ~ scale_y_continuous(limits = c(-31.3883,31.3883)),IsoDecoder == "LysCTT" ~ scale_y_continuous(limits = c(-15.9745,15.9745)),IsoDecoder == "LysTTT" ~ scale_y_continuous(limits = c(-19.8382,19.8382)),IsoDecoder == "Met(e)CAT" ~ scale_y_continuous(limits = c(-16.2335,16.2335)),IsoDecoder == "Met(i)CAT" ~ scale_y_continuous(limits = c(-0.3438,0.3438)),IsoDecoder == "PheGAA" ~ scale_y_continuous(limits = c(-14.7289,14.7289)),IsoDecoder == "ProAGG" ~ scale_y_continuous(limits = c(-3.0497,3.0497)),IsoDecoder == "ProCGG" ~ scale_y_continuous(limits = c(-1.9553,1.9553)),IsoDecoder == "ProTGG" ~ scale_y_continuous(limits = c(-3.0337,3.0337)),IsoDecoder == "SerAGA" ~ scale_y_continuous(limits = c(-10.0222,10.0222)),IsoDecoder == "SerCGA" ~ scale_y_continuous(limits = c(-2.7551,2.7551)),IsoDecoder == "SerGCT" ~ scale_y_continuous(limits = c(-4.1307,4.1307)),IsoDecoder == "SerTGA" ~ scale_y_continuous(limits = c(-2.5718,2.5718)),IsoDecoder == "ThrAGT" ~ scale_y_continuous(limits = c(-4.461,4.461)),IsoDecoder == "ThrCGT" ~ scale_y_continuous(limits = c(-3.6271,3.6271)),IsoDecoder == "ThrTGT" ~ scale_y_continuous(limits = c(-8.8457,8.8457)),IsoDecoder == "TrpCCA" ~ scale_y_continuous(limits = c(-111.8872,111.8872)),IsoDecoder == "TyrGTA" ~ scale_y_continuous(limits = c(-12.813,12.813)),IsoDecoder == "ValAAC" ~ scale_y_continuous(limits = c(-9.6495,9.6495)),IsoDecoder == "ValCAC" ~ scale_y_continuous(limits = c(-2.6356,2.6356)),IsoDecoder == "ValTAC" ~ scale_y_continuous(limits = c(-9.9044,9.9044))))



ggsave("plastid_start_plot.pdf", plot=plastid_start_plot, width=6.7, height=4.5)

ggsave("mitochondrial_start_plot.pdf", plot=mitochondrial_start_plot, width=6.7, height=3.1)

ggsave("nuclear_start_plot.pdf", plot=nuclear_start_plot, width=6.7, height=6.8)

ggsave("plastid_end_plot.pdf", plot=plastid_end_plot, width=6.7, height=4.5)

ggsave("mitochondrial_end_plot.pdf", plot=mitochondrial_end_plot, width=6.7, height=3.1)

ggsave("nuclear_end_plot.pdf", plot=nuclear_end_plot, width=6.7, height=6.8)