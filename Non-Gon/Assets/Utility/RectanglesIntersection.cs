using System.Collections.Generic;
using UnityEngine;

public class RectanglesIntersection : MonoBehaviour
{
    // Function to check if two rectangles intersect
    public static bool RectanglesIntersect(List<Vector3> rect1Points, List<Vector3> rect2Points)
    {
        Debug.Log("RectanglesIntersect called");
        return DoLinesIntersect(rect1Points[0], rect1Points[1], rect2Points[0], rect2Points[1]) ||
            DoLinesIntersect(rect1Points[0], rect1Points[1], rect2Points[1], rect2Points[2]) ||
            DoLinesIntersect(rect1Points[0], rect1Points[1], rect2Points[2], rect2Points[3]) ||
            DoLinesIntersect(rect1Points[0], rect1Points[1], rect2Points[3], rect2Points[0]) ||
            DoLinesIntersect(rect1Points[1], rect1Points[2], rect2Points[0], rect2Points[1]) ||
            DoLinesIntersect(rect1Points[1], rect1Points[2], rect2Points[1], rect2Points[2]) ||
            DoLinesIntersect(rect1Points[1], rect1Points[2], rect2Points[2], rect2Points[3]) ||
            DoLinesIntersect(rect1Points[1], rect1Points[2], rect2Points[3], rect2Points[0]) ||
            DoLinesIntersect(rect1Points[2], rect1Points[3], rect2Points[0], rect2Points[1]) ||
            DoLinesIntersect(rect1Points[2], rect1Points[3], rect2Points[1], rect2Points[2]) ||
            DoLinesIntersect(rect1Points[2], rect1Points[3], rect2Points[2], rect2Points[3]) ||
            DoLinesIntersect(rect1Points[2], rect1Points[3], rect2Points[3], rect2Points[0]) ||
            DoLinesIntersect(rect1Points[3], rect1Points[0], rect2Points[0], rect2Points[1]) ||
            DoLinesIntersect(rect1Points[3], rect1Points[0], rect2Points[1], rect2Points[2]) ||
            DoLinesIntersect(rect1Points[3], rect1Points[0], rect2Points[2], rect2Points[3]) ||
            DoLinesIntersect(rect1Points[3], rect1Points[0], rect2Points[3], rect2Points[0]);
    }

    public static bool DoLinesIntersect(Vector3 A1, Vector3 A2, Vector3 B1, Vector3 B2)
    {
        Debug.Log("DoLinesIntersect called");
        Vector3 directionA = A2 - A1;
        Vector3 directionB = B2 - B1;
        Vector3 crossProduct = Vector3.Cross(directionA, directionB);

        Vector3 lineVector;
        // Check if lines are parallel
        if (crossProduct.sqrMagnitude < 0.001f)
        {
            // Lines are parallel, check if they intersect
            lineVector = B1 - A1;
            float dotProduct = Vector3.Dot(lineVector, directionA);

            // Check if lines coincide
            if (Mathf.Abs(dotProduct) < 0.001f)
            {
                // Lines coincide, return a point on one of the lines
                return true;
            }
            else
            {
                // Lines are parallel but not coincident, no intersection
                return false;
            }
        }

        lineVector = A1 - B1;
        float a = Vector3.Dot(directionA, directionA);
        float b = Vector3.Dot(directionA, directionB);
        float c = Vector3.Dot(directionB, directionB);
        float d = Vector3.Dot(directionA, lineVector);
        float e = Vector3.Dot(directionB, lineVector);

        float denominator = a * c - b * b;

        // Calculate the parameters for the lines
        float s = (b * e - c * d) / denominator;
        float t = (a * e - b * d) / denominator;

        Vector3 closestPointLineA = A1 + directionA * s;
        Vector3 closestPointLineB = B1 + directionB * t;

        return closestPointLineA == closestPointLineB;
    }

}