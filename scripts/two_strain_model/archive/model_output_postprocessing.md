---
title: "Two strain model output post-processing"
author: "James Azam"
date: 'November 03, 2021'
output:
  html_document:
    keep_md: yes
    toc: yes
---








#Summaries

# Unmitigated outbreak


| variant_emergence_day| peak_incidence| outbreak_size|
|---------------------:|--------------:|-------------:|
|                   365|        7669938|      39831273|
|                   151|        7670233|      39832085|
|                   136|        7671305|      39835309|
|                   121|        7677042|      39852560|
|                   106|        7708154|      39944490|
|                    91|        7870066|      40381233|
|                    76|        8538599|      41716889|
|                    61|        9987496|      43478903|
|                    46|       11308950|      44551629|
|                    31|       11981552|      44996561|
|                    16|       12247641|      45158301|
|                     1|       12343005|      45214760|


# Vax only


| variant_emergence_day| npi_intensity| vax_coverage| vax_speed|control_type | peak_incidence| outbreak_size|
|---------------------:|-------------:|------------:|---------:|:------------|--------------:|-------------:|
|                   365|             0|          0.8|      3.75|vax_only     |       202.0010|      997.8043|
|                   151|             0|          0.8|      3.75|vax_only     |       202.0010|     1052.3347|
|                   136|             0|          0.8|      3.75|vax_only     |       202.0010|     1052.3501|
|                   121|             0|          0.8|      3.75|vax_only     |       202.0010|     1052.3593|
|                   106|             0|          0.8|      3.75|vax_only     |       202.0010|     1052.3649|
|                    91|             0|          0.8|      3.75|vax_only     |       202.0010|     1053.9281|
|                    76|             0|          0.8|      3.75|vax_only     |       202.0010|     1067.6837|
|                    61|             0|          0.8|      3.75|vax_only     |       216.1237|     1100.7644|
|                    46|             0|          0.8|      3.75|vax_only     |       250.8561|     1174.6308|
|                    31|             0|          0.8|      3.75|vax_only     |       278.2839|     1365.6541|
|                    16|             0|          0.8|      3.75|vax_only     |       387.4760|     2016.6651|
|                     1|             0|          0.8|      3.75|vax_only     |       988.0279|     5384.9455|
|                   365|             0|          0.8|      1.00|vax_only     |      7682.4972|    77607.0515|
|                   151|             0|          0.8|      1.00|vax_only     |      7741.0849|    78531.2410|
|                   136|             0|          0.8|      1.00|vax_only     |      7771.1818|    79095.4800|
|                   121|             0|          0.8|      1.00|vax_only     |      7832.0341|    80201.1874|
|                   106|             0|          0.8|      1.00|vax_only     |      7964.7238|    82563.9201|
|                    91|             0|          0.8|      1.00|vax_only     |      8290.7234|    88124.8438|
|                    76|             0|          0.8|      1.00|vax_only     |      9229.6414|   102681.7681|
|                    61|             0|          0.8|      1.00|vax_only     |     12589.7244|   145402.5086|
|                    46|             0|          0.8|      1.00|vax_only     |     25973.4715|   286058.8279|
|                    31|             0|          0.8|      1.00|vax_only     |     77531.9978|   791075.5955|
|                    16|             0|          0.8|      1.00|vax_only     |    269298.1405|  2528868.0960|
|                     1|             0|          0.8|      1.00|vax_only     |    838111.5831|  6801353.2658|


# Vax + NPIs 


| variant_emergence_day| npi_intensity| vax_coverage| vax_speed|control_type | peak_incidence| outbreak_size|
|---------------------:|-------------:|------------:|---------:|:------------|--------------:|-------------:|
|                   365|           0.1|          0.8|      3.75|vax+npi      |       134.2293|      624.2046|
|                   151|           0.1|          0.8|      3.75|vax+npi      |       134.2293|      668.4712|
|                   136|           0.1|          0.8|      3.75|vax+npi      |       134.2293|      668.4774|
|                   121|           0.1|          0.8|      3.75|vax+npi      |       134.2293|      668.4810|
|                   106|           0.1|          0.8|      3.75|vax+npi      |       134.2293|      668.4830|
|                    91|           0.1|          0.8|      3.75|vax+npi      |       134.2293|      669.7387|
|                    76|           0.1|          0.8|      3.75|vax+npi      |       134.2293|      680.5211|
|                    61|           0.1|          0.8|      3.75|vax+npi      |       145.5617|      705.3650|
|                    46|           0.1|          0.8|      3.75|vax+npi      |       176.6439|      757.4658|
|                    31|           0.1|          0.8|      3.75|vax+npi      |       195.3710|      880.2994|
|                    16|           0.1|          0.8|      3.75|vax+npi      |       255.9502|     1245.5983|
|                     1|           0.1|          0.8|      3.75|vax+npi      |       534.7534|     2797.1196|
|                   365|           0.1|          0.8|      1.00|vax+npi      |      1666.3793|    16780.7037|
|                   151|           0.1|          0.8|      1.00|vax+npi      |      1666.3793|    17283.0863|
|                   136|           0.1|          0.8|      1.00|vax+npi      |      1715.3186|    17512.7644|
|                   121|           0.1|          0.8|      1.00|vax+npi      |      1737.2100|    17913.0799|
|                   106|           0.1|          0.8|      1.00|vax+npi      |      1779.8740|    18662.1787|
|                    91|           0.1|          0.8|      1.00|vax+npi      |      1871.8764|    20181.6099|
|                    76|           0.1|          0.8|      1.00|vax+npi      |      2098.0779|    23555.7514|
|                    61|           0.1|          0.8|      1.00|vax+npi      |      2756.0476|    31842.0166|
|                    46|           0.1|          0.8|      1.00|vax+npi      |      4863.5978|    54559.5737|
|                    31|           0.1|          0.8|      1.00|vax+npi      |     11743.8840|   124542.8784|
|                    16|           0.1|          0.8|      1.00|vax+npi      |     35879.6803|   365488.4603|
|                     1|           0.1|          0.8|      1.00|vax+npi      |    127110.3533|  1245268.8937|
