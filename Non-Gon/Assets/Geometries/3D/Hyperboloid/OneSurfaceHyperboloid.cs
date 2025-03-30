using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class OneSurfaceHyperboloid : Geometry
{
    public override Vector3 point (float theta, float phi, GeometryCreator variables)
    {
        Vector3 point = new Vector3();

        float coshPhi = (Mathf.Exp(phi) + Mathf.Exp(-phi)) / 2f;
        float sinhPhi = (Mathf.Exp(phi) - Mathf.Exp(-phi)) / 2f;
        
        point[0] = Mathf.Sign(Mathf.Cos(theta))*variables.xradius*Mathf.Abs(Mathf.Cos(theta))*Mathf.Abs(coshPhi);
        point[1] = Mathf.Sign(Mathf.Sin(theta))*variables.yradius*Mathf.Abs(Mathf.Sin(theta))*Mathf.Abs(coshPhi);
        point[2] = Mathf.Sign(sinhPhi) *variables.zradius*Mathf.Abs(sinhPhi);

        return point;
    }

    public override Vector3 normal (float theta, float phi, GeometryCreator variables)
    {
        Vector3 normal = new Vector3();

        float coshPhi = (Mathf.Exp(phi) + Mathf.Exp(-phi)) / 2f;
        float sinhPhi = (Mathf.Exp(phi) - Mathf.Exp(-phi)) / 2f;
        
        normal[0] = Mathf.Sign(Mathf.Cos(theta))*(1f/variables.xradius)*Mathf.Abs(Mathf.Cos(theta))*Mathf.Abs(coshPhi);
        normal[1] = Mathf.Sign(Mathf.Sin(theta))*(1f/variables.yradius)*Mathf.Abs(Mathf.Sin(theta))*Mathf.Abs(coshPhi);
        normal[2] = Mathf.Sign(sinhPhi)*(1f/variables.zradius)*Mathf.Abs(sinhPhi);

        return Vector3.Normalize(normal);
    }
}
