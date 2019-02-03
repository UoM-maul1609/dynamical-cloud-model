function [u1,u2,u3]=problem_test3
% serial algorithm 1
sz=242;


A=2.*eye(sz,sz)+... % diagonal
    [zeros(sz,1),[eye(sz-1,sz-1);zeros(1,sz-1)]]+... % upper
    [[zeros(1,sz-1);0.5.*eye(sz-1,sz-1)],zeros(sz,1)];    % lower
b=ones(sz,1);

u1=A\b;

% serial algorithm 2
%b=2.*ones(sz,1);
%a=0.5.*ones(sz-1,1); % lower diagonal
%a(50)=0.2;
%c=ones(sz-1,1); % upper diagonal
%r=ones(sz,1); % result
a=[-2.4066015339155760       0.40660153391557596       -2.4113484167956374       0.41134841679563738       -2.4161831539184870       0.41618315391848704       -2.4211080643854714       0.42110806438547144       -2.4261255483805959       0.42612554838059591       -2.4312380907317754       0.43123809073177544       -2.4364482646600605       0.43644826466006048       -2.4417587357332025       0.44175873573320246       -2.4471722660291171       0.44717226602911708       -2.4526917185286763       0.45269171852867629       -2.4583200617501406       0.45832006175014062       -2.4640603746410639       0.46406037464106387       -2.4699158517448341       0.46991585174483408       -2.4758898086597201       0.47588980865972008       -2.4819856878101820       0.48198568781018203       -2.4882070645511947       0.48820706455119467       -2.4945576536285374       0.49455765362853743       -2.5010413160192302       0.50104131601923019       -2.5076620661785567       0.50766206617855669       -2.5144240797222484       0.51442407972224835       -2.5213317015744394       0.52133170157443942       -2.5283894546147536       0.52838945461475362       -2.5356020488605866       0.53560204886058660       -2.5429743912235230       0.54297439122352298       -2.5505115958822726       0.55051159588227261       -2.5582189953180405       0.55821899531804053       -2.5661021520621121       0.56610215206211212       -2.5741668712099850       0.57416687120998500       -2.5824192137608701       0.58241921376087014       -2.5908655108468723       0.59086551084687233       -2.5995123789217169       0.59951237892171694       -2.6083667317257415       0.60836673172574152       -2.6174358093607717       0.61743580936077169       -2.6267271859764234       0.62672718597642341       -2.6362487930094316       0.63624879300943160       -2.6460089404666212       0.64600894046662116       -2.6560163397472341       0.65601633974723406       -2.6662801281360915       0.66628012813609150       -2.6768098951119255       0.67680989511192546       -2.6876157106299665       0.68761571062996651       -2.6987081555536721       0.69870815555367205       -2.7100983544286459       0.71009835442864588       -2.7217980108116375       0.72179801081163752       -2.7338194453901909       0.73381944539019095       -2.7461756371535255       0.74617563715352553       -2.7588802679032405       0.75888026790324048       -2.7719477704244131       0.77194777042441309       -2.7853933806730655       0.78539338067306552       -2.7992331943760345       0.79923319437603446       -2.8134842284848469       0.81348422848484692       -2.8281644879760304       0.82816448797603037       -2.8432930385486479       0.84329303854864790       -2.8588900858352622       0.85889008583526216       -2.8749770618174688       0.87497706181746882       -2.8915767192220754       0.89157671922207538       -2.9087132347714304       0.90871323477143040       -2.9264123222715019       0.92641232227150194       -2.9447013566493956       0.94470135664939558       -2.9636095101967843       0.96360951019678431       -2.9831679024440305       0.98316790244403052       -3.0034097652825267        1.0034097652825267       -3.0243706251763349        1.0243706251763349       -3.1597045724320125        1.1597045724320125       -3.1847628944196176        1.1847628944196176       -3.2110982170225775        1.2110982170225775       -3.2382585126015861        1.2382585126015861       -3.2663346832007374        1.2663346832007374       -3.2953385491558418        1.2953385491558418       -3.3259326807372096        1.3259326807372096       -3.3566804709330009        1.3566804709330009       -3.3908245397609695        1.3908245397609695       -3.4253911755267481        1.4253911755267481       -3.4607270977630122        1.4607270977630122       -3.4965963415408297        1.4965963415408297       -3.5370927735866147        1.5370927735866147       -3.5794550155717788        1.5794550155717788       -3.6222782397633999        1.6222782397633999       -3.6659140170168225        1.6659140170168225       -3.7141623990577859        1.7141623990577859       -3.7650384967577306        1.7650384967577306       -3.8184297255240169        1.8184297255240169       -3.8726193944168958        1.8726193944168958       -3.9338129641425441        1.9338129641425441       -3.9939161858733643        1.9939161858733643       -4.0605946682216949        2.0605946682216949       -4.1314080826355513        2.1314080826355513       -4.2063777278496257        2.2063777278496257       -4.2818476408863502        2.2818476408863502       -4.3698669135260886        2.3698669135260886       -4.4606355039948893        2.4606355039948893       -4.5581874541040701        2.5581874541040701       -4.6577354852978194        2.6577354852978194       -4.7696672735802466        2.7696672735802466       -4.8899416189035207        2.8899416189035207       -5.0204010274986537        3.0204010274986537       -5.1626441355443493        3.1626441355443493       -5.3181279331947904        3.3181279331947904       -5.4884658296121263        3.4884658296121263       -5.6758825334111274        3.6758825334111274       -5.8830702419867276        3.8830702419867276       -6.1133182432645832        4.1133182432645832       -6.3706881636716357        4.3706881636716357       -6.6602547882182286        4.6602547882182286       -6.9884429805646544        4.9884429805646544       -7.3635086312026070        5.3635086312026070       -7.7962409755755617        5.7962409755755617       -8.3010150326779026        6.3010150326779035       -8.8974162493947659        6.8974162493947659       -9.6128364738557810        7.6128364738557810       -10.486793759635946        8.4867937596359457       -11.578477670987795        9.5784776709877946       -12.980729270661621        10.980729270661621       -14.847918644493983        12.847918644493983       -17.457038975574804        15.457038975574804       -21.359633392706328        19.359633392706328       -27.834012396738665        25.834012396738665       -40.670399632195590        38.670399632195590       -78.298785872363524        76.298785872363524       -2230.1258575151337        2228.1258575151337       -1.0000000000000000       -1.0000000000000000        0.0000000000000000];
b=[1.0000000000000000        2.8132030678311510       -2.8132030678311510        2.8226968335912757       -2.8226968335912757        2.8323663078369741       -2.8323663078369741        2.8422161287709429       -2.8422161287709429        2.8522510967611927       -2.8522510967611927        2.8624761814635500       -2.8624761814635500        2.8728965293201210       -2.8728965293201210        2.8835174714664058       -2.8835174714664058        2.8943445320582333       -2.8943445320582333        2.9053834370573526       -2.9053834370573526        2.9166401235002803       -2.9166401235002803        2.9281207492821286       -2.9281207492821286        2.9398317034896673       -2.9398317034896673        2.9517796173194402       -2.9517796173194402        2.9639713756203632       -2.9639713756203632        2.9764141291023893       -2.9764141291023893        2.9891153072570749       -2.9891153072570749        3.0020826320384604       -3.0020826320384604        3.0153241323571134       -3.0153241323571134        3.0288481594444967       -3.0288481594444967        3.0426634031488788       -3.0426634031488788        3.0567789092295072       -3.0567789092295072        3.0712040977211732       -3.0712040977211732        3.0859487824470460       -3.0859487824470460        3.1010231917645470       -3.1010231917645470        3.1164379906360802       -3.1164379906360802        3.1322043041242251       -3.1322043041242251        3.1483337424199682       -3.1483337424199682        3.1648384275217394       -3.1648384275217394        3.1817310216937447       -3.1817310216937447        3.1990247578434330       -3.1990247578434330        3.2167334634514839       -3.2167334634514839        3.2348716187215425       -3.2348716187215425        3.2534543719528468       -3.2534543719528468        3.2724975860188632       -3.2724975860188632        3.2920178809332405       -3.2920178809332405        3.3120326794944690       -3.3120326794944690        3.3325602562721830       -3.3325602562721830        3.3536197902238500       -3.3536197902238500        3.3752314212599330       -3.3752314212599330        3.3974163111073441       -3.3974163111073441        3.4201967088572918       -3.4201967088572918        3.4435960216232768       -3.4435960216232768        3.4676388907803810       -3.4676388907803810        3.4923512743070511       -3.4923512743070511        3.5177605358064801       -3.5177605358064801        3.5438955408488271       -3.5438955408488271        3.5707867613461310       -3.5707867613461310        3.5984663887520698       -3.5984663887520698        3.6269684569696947       -3.6269684569696947        3.6563289759520599       -3.6563289759520599        3.6865860770972958       -3.6865860770972958        3.7177801716705252       -3.7177801716705252        3.7499541236349376       -3.7499541236349376        3.7831534384441508       -3.7831534384441508        3.8174264695428599       -3.8174264695428599        3.8528246445430039       -3.8528246445430039        3.8894027132987912       -3.8894027132987912        3.9272190203935686       -3.9272190203935686        3.9663358048880610       -3.9663358048880610        4.0068195305650534       -4.0068195305650534        4.0487412503526707       -4.0487412503526707        4.3194091448640259       -4.3194091448640259        4.3695257888392360       -4.3695257888392360        4.4221964340451541       -4.4221964340451541        4.4765170252031714       -4.4765170252031714        4.5326693664014748       -4.5326693664014748        4.5906770983116836       -4.5906770983116836        4.6518653614744183       -4.6518653614744183        4.7133609418660027       -4.7133609418660027        4.7816490795219382       -4.7816490795219382        4.8507823510534962       -4.8507823510534962        4.9214541955260245       -4.9214541955260245        4.9931926830816593       -4.9931926830816593        5.0741855471732311       -5.0741855471732311        5.1589100311435558       -5.1589100311435558        5.2445564795267998       -5.2445564795267998        5.3318280340336441       -5.3318280340336441        5.4283247981155709       -5.4283247981155709        5.5300769935154612       -5.5300769935154612        5.6368594510480339       -5.6368594510480339        5.7452387888337908       -5.7452387888337908        5.8676259282850891       -5.8676259282850891        5.9878323717467277       -5.9878323717467277        6.1211893364433889       -6.1211893364433889        6.2628161652711016       -6.2628161652711016        6.4127554556992523       -6.4127554556992523        6.5636952817727003       -6.5636952817727003        6.7397338270521772       -6.7397338270521772        6.9212710079897786       -6.9212710079897786        7.1163749082081402       -7.1163749082081402        7.3154709705956407       -7.3154709705956407        7.5393345471604931       -7.5393345471604931        7.7798832378070397       -7.7798832378070397        8.0408020549973074       -8.0408020549973074        8.3252882710886986       -8.3252882710886986        8.6362558663895808       -8.6362558663895808        8.9769316592242525       -8.9769316592242525        9.3517650668222529       -9.3517650668222529        9.7661404839734551       -9.7661404839734551        10.226636486529166       -10.226636486529166        10.741376327343271       -10.741376327343271        11.320509576436457       -11.320509576436457        11.976885961129309       -11.976885961129309        12.727017262405212       -12.727017262405212        13.592481951151122       -13.592481951151122        14.602030065355805       -14.602030065355805        15.794832498789532       -15.794832498789532        17.225672947711562       -17.225672947711562        18.973587519271891       -18.973587519271891        21.156955341975589       -21.156955341975589        23.961458541323246       -23.961458541323246        27.695837288987974       -27.695837288987974        32.914077951149608       -32.914077951149608        40.719266785412664       -40.719266785412664        53.668024793477329       -53.668024793477329        79.340799264391194       -79.340799264391194        154.59757174472705       -154.59757174472705        4458.2517150302674       -4458.2517150302674        0.0000000000000000       -0.0000000000000000        1.0000000000000000];
c=[-0.20999999344348907      -0.40660153391557596        2.4066015339155760      -0.41134841679563738        2.4113484167956374      -0.41618315391848704        2.4161831539184870      -0.42110806438547144        2.4211080643854714      -0.42612554838059591        2.4261255483805959      -0.43123809073177544        2.4312380907317754      -0.43644826466006048        2.4364482646600605      -0.44175873573320246        2.4417587357332025      -0.44717226602911708        2.4471722660291171      -0.45269171852867629        2.4526917185286763      -0.45832006175014062        2.4583200617501406      -0.46406037464106387        2.4640603746410639      -0.46991585174483408        2.4699158517448341      -0.47588980865972008        2.4758898086597201      -0.48198568781018203        2.4819856878101820      -0.48820706455119467        2.4882070645511947      -0.49455765362853743        2.4945576536285374      -0.50104131601923019        2.5010413160192302      -0.50766206617855669        2.5076620661785567      -0.51442407972224835        2.5144240797222484      -0.52133170157443942        2.5213317015744394      -0.52838945461475362        2.5283894546147536      -0.53560204886058660        2.5356020488605866      -0.54297439122352298        2.5429743912235230      -0.55051159588227261        2.5505115958822726      -0.55821899531804053        2.5582189953180405      -0.56610215206211212        2.5661021520621121      -0.57416687120998500        2.5741668712099850      -0.58241921376087014        2.5824192137608701      -0.59086551084687233        2.5908655108468723      -0.59951237892171694        2.5995123789217169      -0.60836673172574152        2.6083667317257415      -0.61743580936077169        2.6174358093607717      -0.62672718597642341        2.6267271859764234      -0.63624879300943160        2.6362487930094316      -0.64600894046662116        2.6460089404666212      -0.65601633974723406        2.6560163397472341      -0.66628012813609150        2.6662801281360915      -0.67680989511192546        2.6768098951119255      -0.68761571062996651        2.6876157106299665      -0.69870815555367205        2.6987081555536721      -0.71009835442864588        2.7100983544286459      -0.72179801081163752        2.7217980108116375      -0.73381944539019095        2.7338194453901909      -0.74617563715352553        2.7461756371535255      -0.75888026790324048        2.7588802679032405      -0.77194777042441309        2.7719477704244131      -0.78539338067306552        2.7853933806730655      -0.79923319437603446        2.7992331943760345      -0.81348422848484692        2.8134842284848469      -0.82816448797603037        2.8281644879760304      -0.84329303854864790        2.8432930385486479      -0.85889008583526216        2.8588900858352622      -0.87497706181746882        2.8749770618174688      -0.89157671922207538        2.8915767192220754      -0.90871323477143040        2.9087132347714304      -0.92641232227150194        2.9264123222715019      -0.94470135664939558        2.9447013566493956      -0.96360951019678431        2.9636095101967843      -0.98316790244403052        2.9831679024440305       -1.0034097652825267        3.0034097652825267       -1.0243706251763349        3.0243706251763349       -1.1597045724320125        3.1597045724320125       -1.1847628944196176        3.1847628944196176       -1.2110982170225775        3.2110982170225775       -1.2382585126015861        3.2382585126015861       -1.2663346832007374        3.2663346832007374       -1.2953385491558418        3.2953385491558418       -1.3259326807372096        3.3259326807372096       -1.3566804709330009        3.3566804709330009       -1.3908245397609695        3.3908245397609695       -1.4253911755267481        3.4253911755267481       -1.4607270977630122        3.4607270977630122       -1.4965963415408297        3.4965963415408297       -1.5370927735866147        3.5370927735866147       -1.5794550155717788        3.5794550155717788       -1.6222782397633999        3.6222782397633999       -1.6659140170168225        3.6659140170168225       -1.7141623990577859        3.7141623990577859       -1.7650384967577306        3.7650384967577306       -1.8184297255240169        3.8184297255240169       -1.8726193944168958        3.8726193944168958       -1.9338129641425441        3.9338129641425441       -1.9939161858733643        3.9939161858733643       -2.0605946682216949        4.0605946682216949       -2.1314080826355513        4.1314080826355513       -2.2063777278496257        4.2063777278496257       -2.2818476408863502        4.2818476408863502       -2.3698669135260886        4.3698669135260886       -2.4606355039948893        4.4606355039948893       -2.5581874541040701        4.5581874541040701       -2.6577354852978194        4.6577354852978194       -2.7696672735802466        4.7696672735802466       -2.8899416189035207        4.8899416189035207       -3.0204010274986537        5.0204010274986537       -3.1626441355443493        5.1626441355443493       -3.3181279331947904        5.3181279331947904       -3.4884658296121263        5.4884658296121263       -3.6758825334111274        5.6758825334111274       -3.8830702419867276        5.8830702419867276       -4.1133182432645832        6.1133182432645832       -4.3706881636716357        6.3706881636716357       -4.6602547882182286        6.6602547882182286       -4.9884429805646544        6.9884429805646544       -5.3635086312026070        7.3635086312026070       -5.7962409755755617        7.7962409755755617       -6.3010150326779035        8.3010150326779026       -6.8974162493947659        8.8974162493947659       -7.6128364738557810        9.6128364738557810       -8.4867937596359457        10.486793759635946       -9.5784776709877946        11.578477670987795       -10.980729270661621        12.980729270661621       -12.847918644493983        14.847918644493983       -15.457038975574804        17.457038975574804       -19.359633392706328        21.359633392706328       -25.834012396738665        27.834012396738665       -38.670399632195590        40.670399632195590       -76.298785872363524        78.298785872363524       -2228.1258575151337        2230.1258575151337        1.0000000000000000        1.0000000000000000];
r=[8.6754109823335927E-131   0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        0.0000000000000000       -0.0000000000000000        3.9525251667299724E-323  -3.9525251667299724E-323   5.3908688375288429E-318  -5.3908688375288429E-318   7.3490830261733734E-313  -7.3490830261733734E-313   9.0255810341636697E-308  -9.0255810341636697E-308   9.9914072518356040E-303  -9.9914072518356040E-303   9.9753665638957889E-298  -9.9753665638957889E-298   8.9871627302536856E-293  -8.9871627302536856E-293   7.3105206674792296E-288  -7.3105206674792296E-288   5.3721235382995148E-283  -5.3721235382995158E-283   3.5682499019930180E-278  -3.5682499019930187E-278   2.1434597338019878E-273  -2.1434597338019874E-273   1.1651048376592365E-268  -1.1651048376592365E-268   5.7338021169912735E-264  -5.7338021169912735E-264   2.5561504680758215E-259  -2.5561504680758212E-259   1.0328416616103270E-254  -1.0328416616103271E-254   3.7846166670006050E-250  -3.7846166670006050E-250   1.2583100488485461E-245  -1.2583100488485461E-245   3.7981092435266075E-241  -3.7981092435266075E-241   1.0413513044707505E-236  -1.0413513044707505E-236   2.5948527429322209E-232  -2.5948527429322209E-232   5.8796168417575000E-228  -5.8796168417575000E-228   1.2121071657328452E-223  -1.2121071657328452E-223   2.2746948080323296E-219  -2.2746948080323296E-219   3.8880284192400628E-215  -3.8880284192400628E-215   6.0561074248510887E-211  -6.0561074248510871E-211   8.6009979025682889E-207  -8.6009979025682889E-207   1.1143667719846881E-202  -1.1143667719846881E-202   1.3178441096765418E-198  -1.3178441096765418E-198   1.4232744751490439E-194  -1.4232744751490436E-194   1.4045390251731962E-190  -1.4045390251731962E-190   1.2671584347600031E-186  -1.2671584347600029E-186   1.0457095058875527E-182  -1.0457095058875529E-182   7.8977782991025289E-179  -7.8977782991025289E-179   5.4618898010365148E-175  -5.4618898010365148E-175   3.4606291864063523E-171  -3.4606291864063523E-171   2.0098841465484939E-167  -2.0098841465484939E-167   1.0705825681396377E-163  -1.0705825681396377E-163   5.2327686334607509E-160  -5.2327686334607509E-160   2.3481933000782666E-156  -2.3481933000782671E-156   9.6795579327012215E-153  -9.6795579327012215E-153   3.6671164098674964E-149  -3.6671164098674964E-149   1.2775244142466470E-145  -1.2775244142466470E-145   3.3364061867491216E-142  -3.3364061867491216E-142   6.7335345808928831E-139  -6.7335345808928831E-139   1.2607984875394352E-135  -1.2607984875394352E-135   2.1908864219329712E-132  -2.1908864219329712E-132   3.5364506059008709E-129  -3.5364506059008709E-129   5.3073253401178842E-126  -5.3073253401178842E-126   7.4066768802980345E-123  -7.4066768802980345E-123   9.6240064751850169E-120  -9.6240064751850169E-120   1.1636073817316988E-116  -1.1636073817316988E-116   1.3082629239664374E-113  -1.3082629239664374E-113   1.3709483735445111E-110  -1.3709483735445114E-110   1.3418345658206037E-107  -1.3418345658206037E-107   1.2251993866187872E-104  -1.2251993866187873E-104   1.0413839591754801E-101  -1.0413839591754801E-101   8.2544423054775438E-099  -8.2544423054775438E-099   6.1178635245235761E-096  -6.1178635245235761E-096   4.2390555005075816E-093  -4.2390555005075816E-093   2.7420351885826877E-090  -2.7420351885826881E-090   1.6567192838464064E-087  -1.6567192838464064E-087   9.3659650480070067E-085  -9.3659650480070067E-085   4.9539368395736627E-082  -4.9539368395736627E-082   2.4517434339249333E-079  -2.4517434339249333E-079   1.1376082742342471E-076  -1.1376082742342471E-076   4.9425553499988228E-074  -4.9425553499988228E-074   2.0119349686747279E-071  -2.0119349686747279E-071   7.6883822055805332E-069  -7.6883822055805332E-069   2.7595743606429110E-066  -2.7595743606429110E-066   9.2832989492621717E-064  -9.2832989492621717E-064   2.9328617953550737E-061  -2.9328617953550740E-061   8.7136949330087025E-059  -8.7136949330087025E-059   2.4388915160538125E-056  -2.4388915160538125E-056   6.4227845736430759E-054  -6.4227845736430759E-054   1.5926201482355513E-051  -1.5926201482355513E-051   3.7200712259740475E-049  -3.7200712259740475E-049   8.1887117929559238E-047  -8.1887117929559228E-047   1.6995049872175997E-044  -1.6995049872175997E-044   3.3276334587121501E-042  -3.3276334587121507E-042   6.1507948344984361E-040  -6.1507948344984361E-040   1.0739944277840552E-037  -1.0739944277840552E-037   1.7727995783799043E-035  -1.7727995783799043E-035   2.7684724218438329E-033  -2.7684724218438323E-033   4.0936487507549798E-031  -4.0936487507549798E-031   5.7368960259461212E-029  -5.7368960259461212E-029   7.6277559981275424E-027  -7.6277559981275424E-027   9.6337439925749109E-025  -9.6337439925749109E-025   1.1574131973796886E-022  -1.1574131973796886E-022   1.3250119215168542E-020  -1.3250119215168542E-020   1.4484788787032272E-018  -1.4484788787032272E-018   1.5161898860242535E-016  -1.5161898860242535E-016   1.5252427288980821E-014  -1.5252427288980821E-014   1.4822496921809338E-012  -1.4822496921809338E-012   1.4024687013229133E-010  -1.4024687013229133E-010   1.3085983429326104E-008  -1.3085983429326102E-008   1.2327453768408796E-006  -1.2327453768408796E-006   1.2348876306927348E-004  -1.2348876306927348E-004   1.5462969259510044E-002  -1.5462969259510044E-002   27.187175382899198       -27.187175382899195       -0.0000000000000000        0.0000000000000000        4.6304736016665446E-002];



u2=tridag_ser(b,a,c,r,sz);


% parallel algorithm https://web.alcf.anl.gov/~zippy/publications/partrid/partrid.html
p=2;
assert(sz./p==round(sz./p));
u3=tridag_par(b,a,c,r,sz,p);


function u=tridag_ser(b,a,c,r,sz)

bet=b(1);
u=zeros(sz,1);
gam=zeros(sz,1);
u(1)=r(1)./bet;
if(bet==0) return;end
for j=2:sz
    gam(j)=c(j-1)./bet;
    bet=b(j)-a(j-1).*gam(j);
    if(bet==0) return;end;
    u(j)=(r(j)-a(j-1).*u(j-1))./bet;
end

for j=sz-1:-1:1
    u(j)=u(j)-gam(j+1).*u(j+1);
end

function u=tridag_par(bn,an,cn,rn,sz,p)

% make sub-arrays+++++++++++++++++++++++++
m=sz./p;
a=zeros(m,p);
b=zeros(m,p);
c=zeros(m,p);
r=zeros(m,p);
for i=1:p
    b(:,i)=bn(1+(i-1).*m:i.*m);
    if(i==1) % lower diagonal
        a(1,i)=0;
        a(2:m,i)=an(1:m-1);
    else
        a(:,i)=an([(i-1).*m:i.*m-1]);
    end
    
    if(i==p) % upper diagonal
        c(1:m-1,i)=cn([1+(i-1)*m:i*m-1]);
        c(m,i)=0.;
    else
        c(:,i)=cn(1+(i-1).*m:i.*m);        
    end
    
    r(:,i)=rn(1+(i-1).*m:i.*m);
end
%------------------------------------------



% part 1: forward elim, back sub and forward sub
xuh=zeros(m,p);
xlh=zeros(m,p);
xr=zeros(m,p);
for j=1:p
    % forward elimination
    xuh(1,j)=c(1,j)/b(1,j);
    xlh(1,j)=r(1,j)/b(1,j);
    for i=2:m
        denom=b(i,j)-a(i,j)*xuh(i-1,j); % mistake had a as c here
        assert(denom~=0);
        xuh(i,j)=c(i,j)/denom;
        xlh(i,j)=(r(i,j)-a(i,j)*xlh(i-1,j))/denom;
    end
    % back substitution
    xr(m,j)=xlh(m,j);
    xlh(m,j)=-xuh(m,j);
    xuh(m,j)=a(m,j)/b(m,j);
    for i=m-1:-1:1
        xr(i,j)=xlh(i,j)-xuh(i,j)*xr(i+1,j);
        xlh(i,j)=-xuh(i,j)*xlh(i+1,j);
        denom=b(i,j)-c(i,j)*xuh(i+1,j);
        denom=sign(denom+1e-60).*max(abs(denom),1e-15);
%         disp(denom);
        assert(denom ~=0);
        xuh(i,j)=a(i,j)/denom; % minus sign wrong in paper
    end
    
    % forward substitution
    xuh(1,j)=-xuh(1,j);
    for i=2:m
        xuh(i,j)=-xuh(i,j)*xuh(i-1,j);
    end
    
end


% part 2: send information back
log2P=log(p)./log(2);
OutData=zeros(8*2^log2P,p);
OutData(1,:)=-1;
OutData(2,:)=xuh(1,:);
OutData(3,:)=xlh(1,:);
OutData(4,:)=-xr(1,:);
OutData(5,:)=xuh(m,:);
OutData(6,:)=xlh(m,:);
OutData(7,:)=-1;
OutData(8,:)=-xr(m,:);

% OutDataLocal=reshape(OutData,[8*p 1]);
% OutDataLocal=zeros(8*2^log2P,1);
% OutData
% for j=1:p
%     for i=0:log2P-1
%         nxfer=8*(2^i);
%         ToProc=1 + mod(j+1-2^i+2*p,p);
%         FromProc=1 + mod(j+2^i,p);
%         [ToProc, FromProc]
%         OutData([1:8]+(FromProc-1)*8,ToProc)=OutData(1:8,FromProc);
%     end
% end
% OutData
d=OutData(1:8,:);

OutData=repmat(d(:),[1 p]);
if(p==1)
    u=xr;
    return;
end

% put outdata into reduced tridiagonal form:
nsig=8*p;
ifirst=8*(p-[p:-1:1])+5
sz2=2*p-2;
x=zeros(m,p);
for j=1:p
    for i=1:2*p-2
        ibase=mod(ifirst(1)+4*(i-1),nsig)
        reduca(i,j)=OutData(ibase,1);
        reducb(i,j)=OutData(ibase+1,1);
        reducc(i,j)=OutData(ibase+2,1);
        reducr(i,j)=OutData(ibase+3,1);
    end
    % solve reduced system
    coeffs(:,j)=tridag_ser(reducb(:,j),reduca(2:end,j),...
        reducc(1:end-1,j),reducr(:,j),sz2);
end


% pick out the appropriate elements of coeffs
for j=1:p
    if(j ~= 1)
        uhcoeff=coeffs(2*j-2,j);
    else
        uhcoeff=0;
    end
    
    if(j ~= p)
        lhcoeff=coeffs(2*j-1,j);
    else
        lhcoeff=0;
    end
    % compute the final solution
    for i=1:m
        x(i,j)=xr(i,j)+uhcoeff*xuh(i,j)+lhcoeff*xlh(i,j);
    end
end



u=reshape(x,[m*p 1]);





