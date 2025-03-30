using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class GeometryCreator : MonoBehaviour
{
    public enum geo_type {_2D, _3D}
    public enum _2Dgeo {Ellipse, Superellipse, Line, Circle, Convex_Line, Convex_Circle, Point};
    public enum _3Dgeo {Ellipsoid, Superellipsoid, Plane, Convex, Point, Cylinder, EllipticParaboloid, 
        OneSurfaceHyperboloid, TwoSurfaceHyperboloid};
    public geo_type Dimension;
    public _2Dgeo Primitive_2D;
    public _3Dgeo Primitive_3D;
    public float xradius;
    public float yradius;
    public float zradius;
    public float e;
    public float radius;
    public float e1;
    public float e2;
    public bool render = false;
    LineRenderer line;
    int segments = 1000;
    public Vector2[] points;
    public Material material_2D;
    public Material material_3D;
    private float prevXRadius;
    private float prevYRadius;
    private float prevRadius;
    private _2Dgeo prevPrimitive2D;



    // Start is called before the first frame update
    void Start()
    {
    
    }

    // Update is called once per frame
    void Update()
    {
        line = gameObject.GetComponent<LineRenderer>();
        if (render)
        {
            if (Dimension == geo_type._2D)
            {
                if (!line)
                {
                    line = gameObject.AddComponent<LineRenderer>();
                    line.SetWidth(4, 4);
                    line.material = material_2D;
                }
                line.useWorldSpace = false;
                if (Primitive_2D == _2Dgeo.Line)
                {
                    line.SetVertexCount(2);
                    line.SetPosition(0, new Vector2(1000f, 0f));
                    line.SetPosition(1, new Vector2(-1000f, 0f));
                } else if (Primitive_2D != prevPrimitive2D || xradius != prevXRadius || yradius != prevYRadius || radius != prevRadius)
                {
                    Create2DGeometry();
                }
            }
        }
        else
        {
            if (line)
            {
                DestroyImmediate(line);
            }
        }
    }
    public void Create2DGeometry()
    {
        line.SetVertexCount(segments + 1);
        float angle = -Mathf.PI / 2f;
        Vector2 point = new Vector2();

        for (int i = 0; i < (segments) + 1; i++)
        {
            switch (Primitive_2D)
            {
                case _2Dgeo.Ellipse:
                    point = Ellipse.point(xradius, yradius, angle);
                    break;
                case _2Dgeo.Superellipse:
                    point = Superellipse.point(xradius, yradius, e, angle);
                    break;

                case _2Dgeo.Circle:
                    point = Ellipse.point(radius, radius, angle);
                    break;
                case _2Dgeo.Convex_Circle:
                    point = Convex_circle.point(angle, radius);
                    break;
                case _2Dgeo.Point:
                    point = Ellipse.point(0.5f, 0.5f, angle);
                    break;
            }
            line.SetPosition(i, point);
            angle += (2 * Mathf.PI / (segments));
        }

        // Update previous values
        prevPrimitive2D = Primitive_2D;
        prevXRadius = xradius;
        prevYRadius = yradius;
        prevRadius = radius;
    }
}


#if UNITY_EDITOR
[CustomEditor(typeof(GeometryCreator)), CanEditMultipleObjects]
public class GeometryCreatorEditor : Editor {
	
	public SerializedProperty 
		Dimension_Prop,
		Primitive_2D_Prop,
		Primitive_3D_Prop,
        xradius_Prop,
        yradius_Prop,
        zradius_Prop,
        e_Prop,
        radius_Prop,
        e1_Prop,
        e2_Prop,
        offset_2D_Prop,
        offset_3D_Prop,
        angle_2D_Prop,
        angle_3D_Prop,
        render_Prop,
        material_2D_Prop,
        material_3D_Prop;


	
	void OnEnable () {
        // Setup the SerializedProperties
        Dimension_Prop = serializedObject.FindProperty ("Dimension");
		Primitive_2D_Prop = serializedObject.FindProperty("Primitive_2D");
		Primitive_3D_Prop = serializedObject.FindProperty ("Primitive_3D");	
        xradius_Prop = serializedObject.FindProperty ("xradius");
        yradius_Prop = serializedObject.FindProperty ("yradius");
        zradius_Prop = serializedObject.FindProperty ("zradius");
        e_Prop = serializedObject.FindProperty ("e");
        radius_Prop = serializedObject.FindProperty ("radius");
        e1_Prop = serializedObject.FindProperty ("e1");
        e2_Prop = serializedObject.FindProperty ("e2");
        offset_2D_Prop = serializedObject.FindProperty ("offset_2D");
        offset_3D_Prop = serializedObject.FindProperty ("offset_3D");
        angle_2D_Prop = serializedObject.FindProperty ("angle_2D");
        angle_3D_Prop = serializedObject.FindProperty ("angle_3D");
        render_Prop = serializedObject.FindProperty ("render");
        material_2D_Prop = serializedObject.FindProperty ("material_2D");
        material_3D_Prop = serializedObject.FindProperty ("material_3D");
    }
	
    public void OnSceneGUI()
    {
        GeometryCreator t = target as GeometryCreator;
        int segments = 1000;
        Vector3[] points = new Vector3[segments+1];
        float angle = -Mathf.PI/2f;
                Vector2 point = new Vector2();
                for (int i = 0; i < (segments)+1; i++)
                {
                    switch(t.Primitive_2D)
                    {
                        case GeometryCreator._2Dgeo.Ellipse:
                            point = Ellipse.point(t.xradius, t.yradius, angle);
                            break;
                        case GeometryCreator._2Dgeo.Superellipse: 
                            point = Superellipse.point(t.xradius, t.yradius, t.e, angle);
                            break;
                    }
                    points[i] = t.transform.TransformPoint(point);
                    angle += (2*Mathf.PI / (segments+1));
                }
        // display an orange disc where the object is
        var color = new Color(1, 0.8f, 0.4f, 1);
        Handles.color = color;
        Handles.DrawPolyLine(points);
        // display object "vlue" in scene
    }
    

	public override void OnInspectorGUI() {
		serializedObject.Update ();
		
        EditorGUILayout.PropertyField( render_Prop, new GUIContent("Render") );
		EditorGUILayout.PropertyField( Dimension_Prop );
        
		
		GeometryCreator.geo_type gt = (GeometryCreator.geo_type)Dimension_Prop.enumValueIndex;
		


		if( (GeometryCreator.geo_type)Dimension_Prop.enumValueIndex == GeometryCreator.geo_type._2D)
        {	
            EditorGUILayout.PropertyField( material_2D_Prop, new GUIContent("Material") );
			EditorGUILayout.PropertyField( Primitive_2D_Prop, new GUIContent("Primitive") );
	
            if( (GeometryCreator._2Dgeo)Primitive_2D_Prop.enumValueIndex == GeometryCreator._2Dgeo.Ellipse)
            {
                EditorGUILayout.PropertyField( xradius_Prop, new GUIContent("X Radius") );
                EditorGUILayout.PropertyField( yradius_Prop, new GUIContent("Y Radius") );
            }
            if( (GeometryCreator._2Dgeo)Primitive_2D_Prop.enumValueIndex == GeometryCreator._2Dgeo.Superellipse)
            {
                EditorGUILayout.PropertyField( xradius_Prop, new GUIContent("X Radius") );
                EditorGUILayout.PropertyField( yradius_Prop, new GUIContent("Y Radius") );
                EditorGUILayout.PropertyField( e_Prop, new GUIContent("e") );
            }
            if( (GeometryCreator._2Dgeo)Primitive_2D_Prop.enumValueIndex == GeometryCreator._2Dgeo.Circle)
            {
                EditorGUILayout.PropertyField( radius_Prop, new GUIContent("Radius") );
            }
            if( (GeometryCreator._2Dgeo)Primitive_2D_Prop.enumValueIndex == GeometryCreator._2Dgeo.Convex_Circle)
            {
                EditorGUILayout.PropertyField( radius_Prop, new GUIContent("Radius") );
            }
        }

		if( (GeometryCreator.geo_type)Dimension_Prop.enumValueIndex == GeometryCreator.geo_type._3D)
        {	
            EditorGUILayout.PropertyField( material_3D_Prop, new GUIContent("Material") );
			EditorGUILayout.PropertyField( Primitive_3D_Prop, new GUIContent("Primitive") );	
            if( (GeometryCreator._3Dgeo)Primitive_3D_Prop.enumValueIndex == GeometryCreator._3Dgeo.Ellipsoid)
            {
                EditorGUILayout.PropertyField( xradius_Prop, new GUIContent("X Radius") );
                EditorGUILayout.PropertyField( yradius_Prop, new GUIContent("Y Radius") );
                EditorGUILayout.PropertyField( zradius_Prop, new GUIContent("Z Radius") );
            }
            if( (GeometryCreator._3Dgeo)Primitive_3D_Prop.enumValueIndex == GeometryCreator._3Dgeo.Superellipsoid)
            {
                EditorGUILayout.PropertyField( xradius_Prop, new GUIContent("X Radius") );
                EditorGUILayout.PropertyField( yradius_Prop, new GUIContent("Y Radius") );
                EditorGUILayout.PropertyField( zradius_Prop, new GUIContent("Z Radius") );
                EditorGUILayout.PropertyField( e1_Prop, new GUIContent("e1") );
                EditorGUILayout.PropertyField( e2_Prop, new GUIContent("e2") );
            }
            if ((GeometryCreator._3Dgeo)Primitive_3D_Prop.enumValueIndex == GeometryCreator._3Dgeo.Cylinder)
            {
                EditorGUILayout.PropertyField(xradius_Prop, new GUIContent("X Radius"));
                EditorGUILayout.PropertyField(yradius_Prop, new GUIContent("Y Radius"));
            }
            if ((GeometryCreator._3Dgeo)Primitive_3D_Prop.enumValueIndex == GeometryCreator._3Dgeo.EllipticParaboloid)
            {
                EditorGUILayout.PropertyField(xradius_Prop, new GUIContent("X Radius"));
                EditorGUILayout.PropertyField(yradius_Prop, new GUIContent("Y Radius"));
            }
            if ((GeometryCreator._3Dgeo)Primitive_3D_Prop.enumValueIndex == GeometryCreator._3Dgeo.OneSurfaceHyperboloid)
            {
                EditorGUILayout.PropertyField(xradius_Prop, new GUIContent("X Radius"));
                EditorGUILayout.PropertyField(yradius_Prop, new GUIContent("Y Radius"));
                EditorGUILayout.PropertyField(zradius_Prop, new GUIContent("Z Radius"));
            }
            if ((GeometryCreator._3Dgeo)Primitive_3D_Prop.enumValueIndex == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid)
            {
                EditorGUILayout.PropertyField(xradius_Prop, new GUIContent("X Radius"));
                EditorGUILayout.PropertyField(yradius_Prop, new GUIContent("Y Radius"));
                EditorGUILayout.PropertyField(zradius_Prop, new GUIContent("Z Radius"));
            }

        }
		
		
		serializedObject.ApplyModifiedProperties ();
	}
}
#endif