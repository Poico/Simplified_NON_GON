using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static Minimum_distance_2D;
using static Minimum_distance_3D;
using static ProximityQuery2D;
using static ProximityQuery3D;
using static Utility;

public class SceneUpdate : MonoBehaviour
{
    public bool isMinimumDistance;
    public bool is2D;
    public bool collision;
    public GameObject obj1;
    public GameObject obj2;
    public GameObject connect; //distance visualizer
    
    LineRenderer line;
    
    // Start is called before the first frame update
    void Start()
    {
        if (connect != null)
        {
            line = connect.GetComponent<LineRenderer>();
            line.SetVertexCount(2);
        }
    }

    // Update is called once per frame
    void Update()
    {
        GeometryCreator col1 = obj1.GetComponent<GeometryCreator>();
        GeometryCreator col2 = obj2.GetComponent<GeometryCreator>();
        if (is2D)
        {
            if (isMinimumDistance)
            {
                Vector2 L = new Vector2();
                Vector2 P = new Vector2();
                if (col1.Primitive_2D == GeometryCreator._2Dgeo.Point & col2.Primitive_2D == GeometryCreator._2Dgeo.Point)
                {
                    L = obj1.transform.position;
                    P = obj2.transform.position;
                }
                if (col1.Primitive_2D == GeometryCreator._2Dgeo.Point & col2.Primitive_2D == GeometryCreator._2Dgeo.Ellipse)
                {
                    Vector2 point = obj1.transform.position;
                    Vector2[] T = Minimum_distance_2D.Point_Ellipse(point, obj2);
                    L = T[0];
                    P = T[1];
                }

                if(col1.Primitive_2D == GeometryCreator._2Dgeo.Ellipse & col2.Primitive_2D == GeometryCreator._2Dgeo.Point)
                {
                    Vector2 point = obj2.transform.position;
                    Vector2[] T = Minimum_distance_2D.Point_Ellipse(point, obj1);
                    L = T[0];
                    P = T[1];
                }

                if(col1.Primitive_2D == GeometryCreator._2Dgeo.Ellipse & col2.Primitive_2D == GeometryCreator._2Dgeo.Ellipse)
                {
                    Vector2[] T = Minimum_distance_2D.Ellipse_ellipse(obj1, obj2);
                    L = T[0];
                    P = T[1];
                    float CA = Minimum_distance_2D.Ellipse_ellipse_distance_of_closest_approach(obj1, obj2);
                }

                if(col1.Primitive_2D == GeometryCreator._2Dgeo.Line & col2.Primitive_2D == GeometryCreator._2Dgeo.Superellipse)
                {
                    Vector2[] LP = Minimum_distance_2D.Superellipse_line(obj1, obj2);
                    L = LP[0];
                    P = LP[1];
                }

                if(col1.Primitive_2D == GeometryCreator._2Dgeo.Superellipse & col2.Primitive_2D == GeometryCreator._2Dgeo.Line)
                {
                    Vector2[] LP = Minimum_distance_2D.Superellipse_line(obj2, obj1);
                    L = LP[0];
                    P = LP[1];
                }





                if(col1.Primitive_2D == GeometryCreator._2Dgeo.Convex_Circle & col2.Primitive_2D == GeometryCreator._2Dgeo.Circle)
                {
                    Vector2[] LP = Minimum_distance_2D.Convex_Circle(obj1, obj2);
                    L = LP[0];
                    P = LP[1];
                }

                if(col1.Primitive_2D == GeometryCreator._2Dgeo.Circle & col2.Primitive_2D == GeometryCreator._2Dgeo.Convex_Circle)
                {
                    Vector2[] LP = Minimum_distance_2D.Convex_Circle(obj2, obj1);
                    L = LP[0];
                    P = LP[1];
                }

                Vector3 Point_ = obj2.transform.InverseTransformPoint(L);
                Vector3 p_ = obj2.transform.InverseTransformPoint(P);

                if(Point_.sqrMagnitude > p_.sqrMagnitude)
                {
                    line.startColor = Color.green;
                    line.endColor = Color.green;
                }
                if(Point_.sqrMagnitude < p_.sqrMagnitude)
                {
                    line.startColor = Color.red;
                    line.endColor = Color.red;
                }
                
                line.SetPosition (0, L);
                line.SetPosition (1, P);

            }
            else
            {
                

            }
        }
        else
        {
            if (isMinimumDistance)
            {
                Vector3 L3 = new Vector3();
                Vector3 P3 = new Vector3();

                if (col1.Primitive_3D == GeometryCreator._3Dgeo.Point & col2.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid)
                {
                    Vector3 point = obj1.transform.position;
                    Vector3[] LP = Minimum_distance_3D.point_Ellipsoid(point, obj2);
                    L3 = LP[0];
                    P3 = LP[1];
                }

                if(col1.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid & col2.Primitive_3D == GeometryCreator._3Dgeo.Point)
                {
                    Vector3 point = obj2.transform.position;
                    Vector3[] LP = Minimum_distance_3D.point_Ellipsoid(point, obj1);
                    L3 = LP[0];
                    P3 = LP[1];
                }





                Vector3 Point_ = obj2.transform.InverseTransformPoint(L3);
                Vector3 p_ = obj2.transform.InverseTransformPoint(P3);

                if(Point_.sqrMagnitude > p_.sqrMagnitude)
                {
                    line.startColor = Color.green;
                    line.endColor = Color.green;
                }
                if(Point_.sqrMagnitude < p_.sqrMagnitude)
                {
                    line.startColor = Color.red;
                    line.endColor = Color.red;
                }
                line.SetPosition (0, L3);
                line.SetPosition (1, P3);
            }
            else
            {
                if (col1.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid & col2.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid)
                {
                    collision = ProximityQuery3D.Ellipsoid_Ellipsoid_Caravantes(obj1, obj2);
                }
                if (col1.Primitive_3D == GeometryCreator._3Dgeo.Cylinder & col2.Primitive_3D == GeometryCreator._3Dgeo.Cylinder)
                {
                    collision = ProximityQuery3D.Cylinder_Cylinder_Chittawadigi(obj1, obj2);
                }
                if ((col1.Primitive_3D == GeometryCreator._3Dgeo.OneSurfaceHyperboloid & col2.Primitive_3D == GeometryCreator._3Dgeo.Plane) ||
                    (col1.Primitive_3D == GeometryCreator._3Dgeo.Plane & col2.Primitive_3D == GeometryCreator._3Dgeo.OneSurfaceHyperboloid) ||
                    (col1.Primitive_3D == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid & col2.Primitive_3D == GeometryCreator._3Dgeo.Plane) ||
                    (col1.Primitive_3D == GeometryCreator._3Dgeo.Plane & col2.Primitive_3D == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid))
                {

                    collision = ProximityQuery3D.Hyperboloid_Plane(obj1, obj2);
                }
                if ((col1.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid & col2.Primitive_3D == GeometryCreator._3Dgeo.EllipticParaboloid) ||
                    (col1.Primitive_3D == GeometryCreator._3Dgeo.EllipticParaboloid & col2.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid))
                {
                    collision = ProximityQuery3D.Ellipsoid_EllipticParaboloid_Brozos(obj1, obj2);
                }
            }
        }
    }
}