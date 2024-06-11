#pragma once

#include <slap/slap.h>

sfloat Xref_data[204] = {5.0,
                         7.0,
                         2.0,
                         -1.4,
                         5.160652560060507,
                         6.839320556675474,
                         1.2130512012101469,
                         -1.8135888664905324,
                         5.248301469968943,
                         6.641397890900794,
                         0.5399269969585725,
                         -2.144864449003066,
                         5.273732403226259,
                         6.413960675321315,
                         -0.031308331812253454,
                         -2.403879862586514,
                         5.246585578869316,
                         6.163779035158997,
                         -0.5116281553266366,
                         -2.599752940659822,
                         5.175449379562446,
                         5.896754823933058,
                         -0.9110958308107685,
                         -2.7407312838589633,
                         5.067948735396972,
                         5.618005512612837,
                         -1.238917052498706,
                         -2.834254942545475,
                         4.930828296188796,
                         5.331941938811838,
                         -1.5034917316648082,
                         -2.8870165334744953,
                         4.77003045942473,
                         5.042340179824062,
                         -1.7124650036165288,
                         -2.9050186462810266,
                         4.590768357721344,
                         4.752407825352388,
                         -1.8727770304511915,
                         -2.893628443152452,
                         4.397593939606499,
                         4.464844933391762,
                         -1.9907113318457161,
                         -2.8576293960600743,
                         4.194461301416077,
                         4.181899956543142,
                         -2.0719414319627116,
                         -2.8012701409123415,
                         3.984785446846749,
                         3.9054209266173876,
                         -2.121575659423859,
                         -2.7283104576027437,
                         3.7714966648816612,
                         3.636902183250262,
                         -2.1441999798779015,
                         -2.642064409739765,
                         3.557090727008504,
                         3.3775269278432885,
                         -2.1439187775852435,
                         -2.5454406983997098,
                         3.3436751114182965,
                         3.128205877872067,
                         -2.1243935342189033,
                         -2.440980301024724,
                         3.1330114656941297,
                         2.8896122888186353,
                         -2.088879380264431,
                         -2.3308914800439156,
                         2.9265545208073664,
                         2.6622136019998983,
                         -2.040259517470838,
                         -2.217082256330822,
                         2.7254876684236526,
                         2.4462999666536716,
                         -1.9810775302034394,
                         -2.1011904505937213,
                         2.5307554109292716,
                         2.242009874046007,
                         -1.9135676196841767,
                         -1.9846114015595697,
                         2.343092889527686,
                         2.04935313028511,
                         -1.8396828083475338,
                         -1.8685234736583753,
                         2.163052690498791,
                         1.8682313831465265,
                         -1.761121172230373,
                         -1.7539114691132918,
                         1.9910291234991955,
                         1.6984564066840826,
                         -1.6793501677615361,
                         -1.6415880601355841,
                         1.827280158820952,
                         1.5397663358507838,
                         -1.5956291258033337,
                         -1.5322133565303908,
                         1.6719472030019182,
                         1.3918400318926178,
                         -1.5110299905773443,
                         -1.4263127226329275,
                         1.5250728842524401,
                         1.2543097479965668,
                         -1.4264563844122196,
                         -1.3242929552880915,
                         1.3866170109673324,
                         1.1267722536465432,
                         -1.342661081289936,
                         -1.226456931712379,
                         1.2564708582464916,
                         1.0087985654273817,
                         -1.260261973126881,
                         -1.1330168326708547,
                         1.134469928951346,
                         0.899942421664707,
                         -1.1797566127760313,
                         -1.0441060425826407,
                         1.0204053274611706,
                         0.7997476283335548,
                         -1.1015354170274778,
                         -0.959789824040401,
                         0.9140338760321995,
                         0.7077543941373491,
                         -1.0258936115519421,
                         -0.8800748598837127,
                         0.8150870955597816,
                         0.6235047635690512,
                         -0.9530419978964161,
                         -0.8049177514822448,
                         0.7232791646444537,
                         0.5465472481283972,
                         -0.8831166204101413,
                         -0.7342325573308361,
                         0.6383139632015288,
                         0.47644074768716427,
                         -0.816187408448358,
                         -0.6678974514938218,
                         0.5598912994563423,
                         0.4127578462669362,
                         -0.7522658664553726,
                         -0.60576057691074,
                         0.4877124120514444,
                         0.3550875592147623,
                         -0.6913118816425852,
                         -0.5476451641327373,
                         0.42148483216843396,
                         0.3030376019213867,
                         -0.6332397160176231,
                         -0.4933539817347754,
                         0.3609266840402354,
                         0.25623624381097576,
                         -0.577923246546348,
                         -0.44267318047344323,
                         0.30577049599831035,
                         0.21433380532435387,
                         -0.5252005142921532,
                         -0.39537558925899446,
                         0.2557665882575768,
                         0.1770038500012115,
                         -0.47487764052251763,
                         -0.3512235172038531,
                         0.21068609797935717,
                         0.1439441185202623,
                         -0.42673216504187494,
                         -0.3099711124151313,
                         0.17032369675533646,
                         0.11487724665804554,
                         -0.3805158594385393,
                         -0.2713663248292037,
                         0.13450005050571714,
                         0.08955130455394529,
                         -0.33595706555384697,
                         -0.23515251725280112,
                         0.10306406686189558,
                         0.0677401903970169,
                         -0.29276260732258425,
                         -0.20106976588576672,
                         0.0758949703847101,
                         0.04924390765463322,
                         -0.25061932222112504,
                         -0.16885588896190687,
                         0.05290424142779979,
                         0.03388875121848737,
                         -0.20919525691708116,
                         -0.1382472397610101,
                         0.03403745006377965,
                         0.021527424324411584,
                         -0.16814057036332172,
                         -0.10897929812050552,
                         0.01927601221861064,
                         0.012039104782798255,
                         -0.12708818654005855,
                         -0.08078709271176107,
                         0.008638890974597507,
                         0.005329475909947972,
                         -0.08565423834020407,
                         -0.05340548474524458,
                         0.002184261870919275,
                         0.0013307345511037905,
                         -0.04343834373336055,
                         -0.02656934243163905,
                         1.1156916756909982e-5,
                         1.5857068697935125e-6,
                         -2.3755349886737942e-5,
                         -1.3634453040884537e-5};

sfloat Uref_data[100] = {
    -7.869487987898531,   -4.135888664905324,  -6.731242042515744,
    -3.312755825125335,   -5.712353287708259,  -2.5901541358344797,
    -4.803198235143832,   -1.9587307807330805, -3.9946767548413185,
    -1.4097834319914113,  -3.2782122168793753, -0.9352365868651165,
    -2.6457467916610216,  -0.5276159092902017, -2.0897327195172055,
    -0.18002112806531131, -1.6031202683466281, 0.11390203128574579,
    -1.1793430139452465,  0.3599904709237788,  -0.8123010011699551,
    0.5635925514773296,   -0.4963422746114752, 0.729596833095979,
    -0.22624320454042407, 0.8624604786297864,  0.0028120229265806273,
    0.9662371134005528,   0.19525243366340164, 1.044603973749857,
    0.35514153954472294,  1.1008882098080828,  0.48619862793592855,
    1.1380922371309359,   0.5918198726739854,  1.1589180573710058,
    0.6750991051926264,   1.1657904903415173,  0.7388481133664283,
    1.1608792790119447,   0.7856163611716078,  1.1461200454508353,
    0.8177100446883694,   1.1232340897770756,  0.8372104195820234,
    1.0937470360519321,   0.8459913522598942,  1.0590063389746323,
    0.8457360616512467,   1.0201976734483607,  0.8379530312228357,
    0.978360235757125,    0.8239910816305508,  0.9344009904152427,
    0.8050536035084973,   0.8891079008821405,  0.7822119574855338,
    0.8431621854223972,   0.7564180547553572,  0.7971496415668827,
    0.7285161365552604,   0.7515710840146792,  0.6992537748627479,
    0.7068519415140865,   0.6692921196178323,  0.6633510583701431,
    0.6392154199298542,   0.6213687458308186,  0.6095398481278744,
    0.5811541277800264,   0.5807216562496212,  0.5429118239796196,
    0.5531646947127509,   0.5068080126133219,  0.5272273225419483,
    0.4729759121444876,   0.5032287376963556,  0.4415207205514137,
    0.481454754806427,    0.4125240478872175,  0.4621630560333569,
    0.3860478758592761,   0.44558793884692294, 0.36213807576402574,
    0.4319445823126274,   0.340827513670344,   0.4214328510145921,
    0.3221387692385986,   0.4142406530404389,  0.3060864920089675,
    0.41054686553759434,  0.29267941640504586, 0.4105238382326318,
    0.28192205408744453,  0.41433948199854476, 0.2738160796651649,
    0.4221589460684351,   0.26836142313605527, 0.43414588383473807,
    0.26555707978598164};