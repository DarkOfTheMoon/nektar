<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D flexible cylinder flow simulation using "MovingBody" module</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_MovBody.xml</parameters>
    <files>
        <file description="Session File">CylFlow_MovBody.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">53.6708</value>
            <value variable="v" tolerance="1e-12">1.60229</value>
            <value variable="w" tolerance="1e-12">0.00915168</value>
            <value variable="p" tolerance="1e-12">165.65</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.89588</value>
            <value variable="v" tolerance="1e-12">0.918691</value>
            <value variable="w" tolerance="1e-12">0.00615171</value>
            <value variable="p" tolerance="1e-12">7.87955</value>
        </metric>
    </metrics>
</test>
