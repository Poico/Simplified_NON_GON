using UnityEngine;

public class CirclePlaneIntersection : MonoBehaviour
{
    // Method to find intersection points between a Plane and a Circle
    public static Vector3[] FindIntersectionPoints(Vector3 planePoint, Vector3 planeNormal, Vector3 circleCenter, Vector3 circleNormal, float radius)
    {
        // Ensure the circle's normal is perpendicular to the plane
        if (Vector3.Dot(planeNormal, circleNormal) >= 0.00001f || Vector3.Dot(planeNormal, circleNormal) <= -0.00001f)
        {
            return null;
        }

        // Project circle center onto the plane
        float distanceToPlane = Vector3.Dot(planeNormal, (circleCenter - planePoint));
        Vector3 projectedCenter = circleCenter - distanceToPlane * planeNormal;

        // Calculate distance from the projected center to the circle center
        float centerDistance = Vector3.Distance(projectedCenter, circleCenter);

        if (centerDistance > radius)
        {
            // No intersection
            return new Vector3[0];
        }
        else if (centerDistance == radius)
        {
            // One intersection point
            return new Vector3[] { projectedCenter };
        }
        else
        {
            // Two intersection points
            float distanceFromProjectedCenter = Mathf.Sqrt(radius * radius - centerDistance * centerDistance);
            Vector3 direction = Vector3.Cross(planeNormal, circleNormal).normalized;

            Vector3 intersectionPoint1 = projectedCenter + direction * distanceFromProjectedCenter;
            Vector3 intersectionPoint2 = projectedCenter - direction * distanceFromProjectedCenter;

            return new Vector3[] { intersectionPoint1, intersectionPoint2 };
        }
    }

    public void TEST()
    {
        float CylinderHeight = 2.0f;

        Vector3 PlanePoint = transform.position;
        Vector3 PlaneNormal = new Vector3(0.76f,0.77f,0.00f);
        Vector3 CircleCenter = transform.position + transform.forward * CylinderHeight;
        Vector3 CircleNormal = transform.forward;
        float CircleRadius = 2f;

        var result = FindIntersectionPoints(PlanePoint, PlaneNormal, CircleCenter, CircleNormal, CircleRadius);

        if (result == null)
        {
            Debug.Log("CIRCLE PLANE INTERSECTION RESULT IS NULL");
            return;
        }

        if (result.Length == 0) 
        {
            Debug.Log("CIRCLE DOES NOT INTERSECT PLANE");
            return;
        }

        Debug.Log("CIRCLE PLANE INTERSECTION RESULTED IN " + result.Length + " POINTS");
        for (int i = 0; i < result.Length; i++)
            Debug.Log("CIRCLE PLANE INTERSECTION POINT #"+i+": " + result[i]);
    }

    private void Start()
    {
        TEST();
    }
}