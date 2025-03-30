using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class EllipticParaboloid : Geometry
{
    public override Vector3 point (float theta, float phi, GeometryCreator variables)
    {
        float x = 1.5f * variables.xradius * Mathf.Cos(theta) * Mathf.Sin(phi);
        float y = 1.5f * variables.yradius * Mathf.Sin(theta) * Mathf.Sin(phi);
        float z = 1.5f * ((x * x) / (variables.xradius * variables.xradius) + (y * y) / (variables.yradius * variables.yradius)); // z = (x^2 / a^2) + (y^2 / b^2)

        return new Vector3(x, y, z);
    }

    public override Vector3 normal(float theta, float phi, GeometryCreator variables)
    {
        float x = 1.5f * variables.xradius * Mathf.Cos(theta) * Mathf.Sin(phi);
        float y = 1.5f * variables.yradius * Mathf.Sin(theta) * Mathf.Sin(phi);
        float z = 1.5f * ((x * x) / (variables.xradius * variables.xradius) + (y * y) / (variables.yradius * variables.yradius));

        // Partial derivatives of the parametric equation with respect to theta and phi
        float dx_dtheta = -variables.xradius * Mathf.Sin(theta) * Mathf.Sin(phi);
        float dy_dtheta = variables.yradius * Mathf.Cos(theta) * Mathf.Sin(phi);
        float dz_dtheta = 2f * x / (variables.xradius * variables.xradius) * dx_dtheta + 2f * y / (variables.yradius * variables.yradius) * dy_dtheta;

        float dx_dphi = variables.xradius * Mathf.Cos(theta) * Mathf.Cos(phi);
        float dy_dphi = variables.yradius * Mathf.Sin(theta) * Mathf.Cos(phi);
        float dz_dphi = 2f * x / (variables.xradius * variables.xradius) * dx_dphi + 2f * y / (variables.yradius * variables.yradius) * dy_dphi;

        // Calculate the cross product of the partial derivatives to get the normal vector
        Vector3 normal = new Vector3(
            dz_dtheta * dy_dphi - dz_dphi * dy_dtheta,
            dz_dphi * dx_dtheta - dz_dtheta * dx_dphi,
            dx_dphi * dy_dtheta - dx_dtheta * dy_dphi
        ).normalized;

        return normal;
    }
}
