function C = cmap_mandrill(N)
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2019 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Rana El Khoury, 2019

cmap = [
    0.0803    0.0095    0.2836
    0.0833    0.0164    0.2942
    0.0864    0.0239    0.3048
    0.0893    0.0318    0.3154
    0.0922    0.0403    0.3262
    0.0949    0.0487    0.3371
    0.0976    0.0566    0.3480
    0.1002    0.0642    0.3591
    0.1027    0.0715    0.3703
    0.1050    0.0786    0.3815
    0.1073    0.0855    0.3929
    0.1094    0.0923    0.4044
    0.1114    0.0990    0.4159
    0.1133    0.1055    0.4276
    0.1151    0.1119    0.4393
    0.1166    0.1183    0.4511
    0.1178    0.1249    0.4626
    0.1185    0.1316    0.4738
    0.1187    0.1386    0.4847
    0.1183    0.1458    0.4952
    0.1174    0.1532    0.5051
    0.1158    0.1609    0.5144
    0.1137    0.1688    0.5231
    0.1111    0.1769    0.5310
    0.1079    0.1852    0.5383
    0.1044    0.1936    0.5449
    0.1005    0.2021    0.5509
    0.0963    0.2106    0.5563
    0.0919    0.2191    0.5612
    0.0873    0.2276    0.5657
    0.0825    0.2361    0.5699
    0.0777    0.2445    0.5737
    0.0727    0.2528    0.5773
    0.0677    0.2610    0.5807
    0.0627    0.2692    0.5839
    0.0576    0.2773    0.5870
    0.0525    0.2853    0.5899
    0.0475    0.2933    0.5928
    0.0424    0.3011    0.5956
    0.0374    0.3089    0.5983
    0.0328    0.3167    0.6010
    0.0286    0.3244    0.6036
    0.0248    0.3320    0.6063
    0.0215    0.3396    0.6089
    0.0186    0.3471    0.6115
    0.0161    0.3546    0.6141
    0.0141    0.3620    0.6166
    0.0131    0.3694    0.6191
    0.0134    0.3768    0.6214
    0.0152    0.3841    0.6236
    0.0185    0.3913    0.6257
    0.0239    0.3985    0.6276
    0.0314    0.4056    0.6294
    0.0413    0.4126    0.6310
    0.0527    0.4196    0.6325
    0.0647    0.4264    0.6338
    0.0774    0.4332    0.6350
    0.0905    0.4399    0.6361
    0.1039    0.4465    0.6370
    0.1176    0.4530    0.6379
    0.1316    0.4593    0.6387
    0.1456    0.4656    0.6395
    0.1596    0.4718    0.6402
    0.1738    0.4779    0.6410
    0.1878    0.4839    0.6418
    0.2017    0.4898    0.6427
    0.2156    0.4957    0.6436
    0.2293    0.5015    0.6447
    0.2428    0.5072    0.6458
    0.2562    0.5129    0.6470
    0.2694    0.5185    0.6484
    0.2824    0.5241    0.6498
    0.2953    0.5297    0.6514
    0.3080    0.5352    0.6530
    0.3205    0.5407    0.6548
    0.3329    0.5462    0.6567
    0.3451    0.5517    0.6587
    0.3572    0.5572    0.6608
    0.3692    0.5627    0.6631
    0.3810    0.5681    0.6654
    0.3927    0.5736    0.6678
    0.4042    0.5791    0.6703
    0.4157    0.5845    0.6730
    0.4271    0.5900    0.6757
    0.4383    0.5955    0.6785
    0.4495    0.6010    0.6814
    0.4606    0.6066    0.6844
    0.4716    0.6121    0.6875
    0.4826    0.6177    0.6906
    0.4934    0.6232    0.6939
    0.5042    0.6288    0.6972
    0.5149    0.6345    0.7007
    0.5256    0.6401    0.7042
    0.5363    0.6458    0.7077
    0.5468    0.6514    0.7113
    0.5573    0.6572    0.7150
    0.5678    0.6629    0.7188
    0.5782    0.6687    0.7226
    0.5885    0.6745    0.7264
    0.5988    0.6804    0.7303
    0.6091    0.6863    0.7343
    0.6193    0.6922    0.7383
    0.6295    0.6982    0.7424
    0.6396    0.7042    0.7466
    0.6497    0.7103    0.7508
    0.6598    0.7164    0.7551
    0.6698    0.7225    0.7594
    0.6799    0.7286    0.7638
    0.6899    0.7349    0.7683
    0.6999    0.7411    0.7728
    0.7098    0.7474    0.7774
    0.7198    0.7537    0.7820
    0.7297    0.7601    0.7867
    0.7396    0.7665    0.7915
    0.7495    0.7730    0.7963
    0.7594    0.7795    0.8012
    0.7693    0.7860    0.8062
    0.7792    0.7926    0.8113
    0.7890    0.7993    0.8164
    0.7989    0.8060    0.8215
    0.8087    0.8128    0.8268
    0.8185    0.8196    0.8321
    0.8283    0.8264    0.8376
    0.8381    0.8333    0.8431
    0.8479    0.8403    0.8486
    0.8576    0.8474    0.8543
    0.8673    0.8545    0.8601
    0.8770    0.8616    0.8661
    0.8819    0.8652    0.8691
    0.8771    0.8560    0.8630
    0.8727    0.8468    0.8570
    0.8686    0.8374    0.8510
    0.8648    0.8281    0.8450
    0.8612    0.8186    0.8391
    0.8577    0.8092    0.8332
    0.8544    0.7997    0.8274
    0.8512    0.7902    0.8216
    0.8481    0.7807    0.8159
    0.8451    0.7712    0.8103
    0.8422    0.7616    0.8047
    0.8394    0.7521    0.7991
    0.8367    0.7425    0.7936
    0.8341    0.7329    0.7882
    0.8316    0.7233    0.7828
    0.8291    0.7137    0.7775
    0.8267    0.7040    0.7723
    0.8243    0.6944    0.7671
    0.8221    0.6847    0.7620
    0.8198    0.6750    0.7569
    0.8177    0.6653    0.7519
    0.8156    0.6555    0.7469
    0.8135    0.6457    0.7420
    0.8115    0.6359    0.7372
    0.8095    0.6261    0.7324
    0.8076    0.6162    0.7276
    0.8057    0.6063    0.7230
    0.8039    0.5963    0.7183
    0.8021    0.5863    0.7138
    0.8003    0.5762    0.7092
    0.7986    0.5661    0.7048
    0.7969    0.5560    0.7004
    0.7952    0.5457    0.6960
    0.7936    0.5354    0.6917
    0.7920    0.5251    0.6874
    0.7904    0.5147    0.6831
    0.7889    0.5042    0.6787
    0.7874    0.4936    0.6743
    0.7859    0.4829    0.6699
    0.7846    0.4721    0.6654
    0.7832    0.4612    0.6608
    0.7819    0.4501    0.6562
    0.7807    0.4389    0.6515
    0.7794    0.4276    0.6466
    0.7782    0.4162    0.6417
    0.7770    0.4046    0.6366
    0.7759    0.3928    0.6313
    0.7747    0.3809    0.6258
    0.7736    0.3688    0.6201
    0.7725    0.3565    0.6142
    0.7713    0.3441    0.6079
    0.7701    0.3314    0.6013
    0.7689    0.3186    0.5943
    0.7675    0.3057    0.5869
    0.7661    0.2926    0.5789
    0.7645    0.2794    0.5704
    0.7628    0.2663    0.5613
    0.7607    0.2532    0.5515
    0.7584    0.2402    0.5411
    0.7557    0.2275    0.5300
    0.7526    0.2153    0.5182
    0.7490    0.2035    0.5059
    0.7450    0.1923    0.4931
    0.7405    0.1817    0.4799
    0.7355    0.1717    0.4665
    0.7302    0.1623    0.4529
    0.7245    0.1535    0.4391
    0.7185    0.1452    0.4253
    0.7122    0.1374    0.4115
    0.7056    0.1300    0.3977
    0.6988    0.1229    0.3840
    0.6918    0.1163    0.3704
    0.6846    0.1099    0.3569
    0.6773    0.1037    0.3435
    0.6699    0.0977    0.3301
    0.6623    0.0919    0.3169
    0.6546    0.0863    0.3038
    0.6469    0.0807    0.2909
    0.6391    0.0751    0.2780
    0.6313    0.0694    0.2653
    0.6234    0.0636    0.2526
    0.6154    0.0577    0.2400
    0.6074    0.0517    0.2276
    0.5994    0.0456    0.2153
    0.5913    0.0395    0.2031
    0.5830    0.0335    0.1910
    0.5748    0.0280    0.1790
    0.5664    0.0230    0.1672
    0.5580    0.0185    0.1554
    0.5494    0.0146    0.1438
    0.5407    0.0112    0.1324
    0.5319    0.0083    0.1211
    0.5230    0.0061    0.1100
    0.5139    0.0044    0.0991
    0.5047    0.0034    0.0883
    0.4952    0.0031    0.0779
    0.4856    0.0034    0.0677
    0.4758    0.0044    0.0578
    0.4657    0.0060    0.0484
    0.4554    0.0082    0.0395
    0.4449    0.0109    0.0318
    0.4342    0.0140    0.0257
    0.4233    0.0172    0.0210
    0.4124    0.0204    0.0174
    0.4013    0.0234    0.0148
    0.3903    0.0261    0.0130
    0.3793    0.0284    0.0118
    0.3683    0.0302    0.0110
    0.3573    0.0316    0.0106
    0.3465    0.0326    0.0104
    0.3356    0.0332    0.0103
    0.3248    0.0335    0.0103
    0.3141    0.0336    0.0102
    0.3034    0.0334    0.0101
    0.2927    0.0330    0.0098
    0.2820    0.0324    0.0095
    0.2714    0.0316    0.0092
    0.2607    0.0306    0.0088
    0.2501    0.0295    0.0084
    0.2396    0.0283    0.0079
    0.2290    0.0269    0.0074
    0.2184    0.0254    0.0068
    0.2079    0.0238    0.0062
    0.1973    0.0222    0.0057
    0.1868    0.0205    0.0051
    0.1762    0.0188    0.0045
    0.1656    0.0170    0.0039];

P = size(cmap,1);

if nargin < 1
   N = P;
end

N = min(N,P);
C = interp1(1:P, cmap, linspace(1,P,N), 'linear');