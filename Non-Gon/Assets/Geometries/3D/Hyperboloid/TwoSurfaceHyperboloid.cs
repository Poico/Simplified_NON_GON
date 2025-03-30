using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public class TwoSurfaceHyperboloid : Geometry
{
    public override Vector3 point(float theta, float phi, GeometryCreator variables)
    {
        Vector3 point = new Vector3();

        float coshPhi = (Mathf.Exp(phi) + Mathf.Exp(-phi)) / 2f;
        float sinhPhi = (Mathf.Exp(phi) - Mathf.Exp(-phi)) / 2f;

        point[0] = variables.xradius * Mathf.Cos(theta) * Mathf.Abs(sinhPhi);
        point[1] = variables.yradius * Mathf.Sin(theta) * Mathf.Abs(sinhPhi);
        point[2] = variables.zradius * coshPhi;

        return point;
    }

    public override Vector3 normal (float theta, float phi, GeometryCreator variables)
    {
        Vector3 normal = new Vector3();

        float coshPhi = (Mathf.Exp(phi) + Mathf.Exp(-phi)) / 2f;
        float sinhPhi = (Mathf.Exp(phi) - Mathf.Exp(-phi)) / 2f;
        
        normal[0] = 1f / variables.xradius * Mathf.Cos(theta) * sinhPhi;
        normal[1] = 1f / variables.yradius * Mathf.Sin(theta) * sinhPhi;
        normal[2] = 1f / variables.zradius * coshPhi * Mathf.Sign(sinhPhi);

        return Vector3.Normalize(normal);
    }
}
