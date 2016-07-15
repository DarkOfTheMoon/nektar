<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Hex Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_Kovasnay_SubStep.xml</parameters>
    <files>
        <file description="Session File">Hex_Kovasnay_SubStep.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.0484179 </value>
            <value variable="v" tolerance="1e-6">0.00261627</value>
            <value variable="w" tolerance="1e-6">0.000820649</value>
            <value variable="p" tolerance="1e-6">0.00903438 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0616737 </value>
            <value variable="v" tolerance="1e-6">0.00743594</value>
            <value variable="w" tolerance="1e-6">0.00570445</value>
            <value variable="p" tolerance="1e-6">0.0443199 </value>
        </metric>
    </metrics>
</test>
