using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class Renderer_3d : MonoBehaviour
{
    
    int segments = 30;
    Mesh mesh;
    float phi;
    float theta;
    private float prevXRadius;
    private float prevYRadius;
    private float prevZRadius;
    private float preve1;
    private float preve2;


    void Start()
    {
    }

    void Update()
    {
        GeometryCreator variables = gameObject.GetComponent<GeometryCreator>();
        if (variables.xradius != prevXRadius || variables.yradius != prevYRadius || variables.zradius != prevZRadius || variables.e1 != preve1 || variables.e2 != preve2)
        {
            Create3DGeometry();
        }
    }

    private void Create3DGeometry()
    {
        GeometryCreator variables = gameObject.GetComponent<GeometryCreator>();
        Geometry geometry = null;
        switch (variables.Primitive_3D)
        {
            case GeometryCreator._3Dgeo.Ellipsoid:
                geometry = new Ellipsoid();
                break;
            case GeometryCreator._3Dgeo.Superellipsoid:
                geometry = new Superellipsoid();
                break;
            case GeometryCreator._3Dgeo.Convex:
                geometry = new Parameters_3d();
                break;
            case GeometryCreator._3Dgeo.Cylinder:
                geometry = new Cylinder();
                break;
            case GeometryCreator._3Dgeo.EllipticParaboloid:
                geometry = new EllipticParaboloid();
                break;
            case GeometryCreator._3Dgeo.OneSurfaceHyperboloid:
                geometry = new OneSurfaceHyperboloid();
                break;
            case GeometryCreator._3Dgeo.TwoSurfaceHyperboloid:
                geometry = new TwoSurfaceHyperboloid();
                break;
        }

        float theta = -Mathf.PI;
        float phi = -Mathf.PI / 2 + Mathf.PI / (segments + 1);
        Vector3 point;
        if (variables.Primitive_3D == GeometryCreator._3Dgeo.Convex)
        {
            point = geometry.point(theta, phi + Mathf.PI / 2f, variables);
        }
        else
        {
            point = geometry.point(theta, phi, variables);
        }

        Vector3 normal = geometry.normal(theta, phi, variables);

        mesh = GetComponent<MeshFilter>().mesh;

        int numberOfVertices = 2 + (segments - 1) * (segments + 1);
        Vector3[] vertices = new Vector3[numberOfVertices];
        Vector3[] normals = new Vector3[numberOfVertices];
        int[] triangles = new int[6 * (segments + 1) * (segments - 1)];

        if (variables.Primitive_3D != GeometryCreator._3Dgeo.Cylinder)
        {
            vertices[0] = geometry.point(-Mathf.PI, -Mathf.PI / 2f, variables);
            vertices[vertices.Length - 1] = geometry.point(Mathf.PI, Mathf.PI / 2f, variables);
            normals[0] = geometry.normal(-Mathf.PI + 0.1f, -Mathf.PI / 2f + 0.1f, variables);
            normals[vertices.Length - 1] = geometry.normal(Mathf.PI - 0.1f, Mathf.PI / 2f - 0.1f, variables);
        }
        for (int i = 0; i < (segments); i++)
        {
            vertices[i + 1] = geometry.point(theta, phi, variables);
            normals[i + 1] = geometry.normal(theta, phi, variables);
            if (variables.Primitive_3D != GeometryCreator._3Dgeo.EllipticParaboloid & variables.Primitive_3D != GeometryCreator._3Dgeo.TwoSurfaceHyperboloid &
                variables.Primitive_3D != GeometryCreator._3Dgeo.OneSurfaceHyperboloid & variables.Primitive_3D != GeometryCreator._3Dgeo.Cylinder)
            {
                triangles[i * 3 + 2] = 0;
                triangles[i * 3 + 1] = i + 1;
                triangles[i * 3] = i + 2;
            }
            if (variables.Primitive_3D == GeometryCreator._3Dgeo.Cylinder)
            {
                triangles[i * 3 + 2] = 31;
                triangles[i * 3 + 1] = i + 1;
                triangles[i * 3] = i + 2;
            }
            theta = theta + 2f * Mathf.PI / (segments + 1);
        }
        if (variables.Primitive_3D != GeometryCreator._3Dgeo.TwoSurfaceHyperboloid & variables.Primitive_3D != GeometryCreator._3Dgeo.OneSurfaceHyperboloid &
            variables.Primitive_3D != GeometryCreator._3Dgeo.Cylinder)
        {
            triangles[(segments) * 3 + 2] = 0;
            triangles[(segments) * 3 + 1] = segments + 1;
            triangles[(segments) * 3] = 1;
        }

        theta = -Mathf.PI;

        for (int i = 1; i < (segments); i++)
        {
            for (int j = 0; j < (segments + 1); j++)
            {
                Vector3 point_ = geometry.point(theta, phi, variables);
                Vector3 normal_ = geometry.normal(theta, phi, variables);

                vertices[(i - 1) * (segments + 1) + j + 1] = point_;
                normals[(i - 1) * (segments + 1) + j + 1] = normal_;

                if (i < segments - 1)
                {
                    if (j < segments)
                    {
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + j * 6] = (i - 1) * (segments + 1) + j + 2;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + j * 6 + 1] = (i) * (segments + 1) + j + 1;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + j * 6 + 2] = (i - 1) * (segments + 1) + j + 1;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + j * 6 + 3] = (i - 1) * (segments + 1) + j + 2;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + j * 6 + 4] = (i) * (segments + 1) + j + 2;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + j * 6 + 5] = (i) * (segments + 1) + j + 1;
                    }
                    else
                    {
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + segments * 6] = (i - 1) * (segments + 1) + 1;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + segments * 6 + 1] = (i) * (segments + 1) + segments + 1;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + segments * 6 + 2] = (i - 1) * (segments + 1) + segments + 1;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + segments * 6 + 3] = (i - 1) * (segments + 1) + 1;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + segments * 6 + 4] = (i) * (segments + 1) + 1;
                        triangles[3 * (segments + 1) + (i - 1) * (segments + 1) * 6 + segments * 6 + 5] = (i) * (segments + 1) + segments + 1;
                    }
                }

                theta += (2f * Mathf.PI / (segments + 1));
            }
            theta = -Mathf.PI;
            phi += (Mathf.PI / (segments + 1));
        }

        // Update previous values
        prevXRadius = variables.xradius;
        prevYRadius = variables.yradius;
        prevZRadius = variables.zradius;
        preve1 = variables.e1;
        preve2 = variables.e2;

        if (variables.Primitive_3D != GeometryCreator._3Dgeo.TwoSurfaceHyperboloid & variables.Primitive_3D != GeometryCreator._3Dgeo.EllipticParaboloid &
                variables.Primitive_3D != GeometryCreator._3Dgeo.OneSurfaceHyperboloid)
        {
            if (variables.Primitive_3D == GeometryCreator._3Dgeo.Cylinder)
            {
                for (int i = 0; i < segments; i++)
                {
                    triangles[triangles.Length - 1 - i * 3] = vertices.Length - i - 2;
                    triangles[triangles.Length - 1 - i * 3 - 1] = vertices.Length - 32;
                    triangles[triangles.Length - 1 - i * 3 - 2] = vertices.Length - i - 1;
                }
            }
            else
            {
                for (int i = 0; i < segments; i++)
                {
                    triangles[triangles.Length - 1 - i * 3] = vertices.Length - i - 2;
                    triangles[triangles.Length - 1 - i * 3 - 1] = vertices.Length - 1;
                    triangles[triangles.Length - 1 - i * 3 - 2] = vertices.Length - i - 1;
                }
            }

            if (variables.Primitive_3D != GeometryCreator._3Dgeo.Cylinder)
            {
                triangles[triangles.Length - 1 - (segments) * 3 - 1] = vertices.Length - segments + 1 - 2;
                triangles[triangles.Length - 1 - (segments) * 3 - 2] = vertices.Length - 2;
            }
        }
        triangles[triangles.Length - 1 - (segments) * 3] = vertices.Length - 1;

        if (variables.Primitive_3D == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid)
        {

            // Create arrays to hold positive and negative vertices, triangles, and normals
            Vector3[] mirroredVertices = new Vector3[vertices.Length * 2];
            Vector3[] mirroredNormals = new Vector3[normals.Length * 2];
            int[] mirroredTriangles = new int[triangles.Length * 2];

            // Copy positive vertices, triangles, and normals to the new arrays
            vertices.CopyTo(mirroredVertices, 0);
            normals.CopyTo(mirroredNormals, 0);
            triangles.CopyTo(mirroredTriangles, 0);

            // Mirror vertices and normals along the z-axis for the negative side
            for (int i = 0; i < vertices.Length; i++)
            {
                mirroredVertices[i + vertices.Length] = new Vector3(vertices[i].x, vertices[i].y, -vertices[i].z);
                mirroredNormals[i + normals.Length] = new Vector3(normals[i].x, normals[i].y, -normals[i].z);
            }

            // Offset the triangles indices for the negative side
            for (int i = 0; i < triangles.Length; i += 3)
            {
                mirroredTriangles[i + 0 + triangles.Length] = triangles[i + 2] + vertices.Length;
                mirroredTriangles[i + 1 + triangles.Length] = triangles[i + 1] + vertices.Length;
                mirroredTriangles[i + 2 + triangles.Length] = triangles[i] + vertices.Length;
            }

            // Assign mirrored vertices, triangles, and normals to the mesh
            mesh.vertices = mirroredVertices;
            mesh.triangles = mirroredTriangles;
            mesh.normals = mirroredNormals;
        }
        else
        {
            mesh.vertices = vertices;
            mesh.triangles = triangles;
            mesh.normals = normals;
        }
    }
}
