using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Superellipsoid : Geometry
{

    public override Vector3 point (float theta, float phi, GeometryCreator variables)
    {
        Vector3 point = new Vector3();

        point[0] = Mathf.Sign(Mathf.Cos(theta))*variables.xradius*Mathf.Pow(Mathf.Abs(Mathf.Cos(theta)),variables.e1)*Mathf.Pow(Mathf.Abs(Mathf.Cos(phi)),variables.e2);
        point[1] = Mathf.Sign(Mathf.Sin(theta))*variables.yradius*Mathf.Pow(Mathf.Abs(Mathf.Sin(theta)),variables.e1)*Mathf.Pow(Mathf.Abs(Mathf.Cos(phi)),variables.e2);
        point[2] = Mathf.Sign(Mathf.Sin(phi)) * variables.zradius*Mathf.Pow(Mathf.Abs(Mathf.Sin(phi)),variables.e2);

        return point;
    }

    public override Vector3 normal (float theta, float phi, GeometryCreator variables)
    {
        Vector3 normal = new Vector3();

        normal[0] = Mathf.Sign(Mathf.Cos(theta))*(1f/variables.xradius)*Mathf.Pow(Mathf.Abs(Mathf.Cos(theta)),2f-variables.e1)*Mathf.Pow(Mathf.Abs(Mathf.Cos(phi)),2f-variables.e2);
        normal[1] = Mathf.Sign(Mathf.Sin(theta))*(1f/variables.yradius)*Mathf.Pow(Mathf.Abs(Mathf.Sin(theta)),2f-variables.e1)*Mathf.Pow(Mathf.Abs(Mathf.Cos(phi)),2f-variables.e2);
        normal[2] = Mathf.Sign(Mathf.Sin(phi))*(1f/variables.zradius)*Mathf.Pow(Mathf.Abs(Mathf.Cos(phi)),2f-variables.e2);

        return Vector3.Normalize(normal);
    }
}
