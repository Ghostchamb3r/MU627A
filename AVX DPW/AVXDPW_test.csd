<CsoundSynthesizer>
<CsOptions>
-odac
--opcode-lib=AVXDPW.dll

</CsOptions>
<CsInstruments>

sr = 44100
ksmps = 64
nchnls = 1
0dbfs = 1


instr 1 
icps cpsmidi 
iamp ampmidi 1


asig AVXDPW icps, iamp



outs asig
endin


</CsInstruments>
<CsScore>
</CsScore>
</CsoundSynthesizer>




















<bsbPanel>
 <label>Widgets</label>
 <objectName/>
 <x>100</x>
 <y>100</y>
 <width>320</width>
 <height>240</height>
 <visible>true</visible>
 <uuid/>
 <bgcolor mode="background">
  <r>240</r>
  <g>240</g>
  <b>240</b>
 </bgcolor>
</bsbPanel>
<bsbPresets>
</bsbPresets>
