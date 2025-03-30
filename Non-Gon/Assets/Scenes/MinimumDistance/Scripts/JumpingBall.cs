using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PendulumMovement : MonoBehaviour
{
    public float targetDistance = 1000f; // Distance to move to the left
    public float arcHeight = 100f; // Height of the arc
    public float duration = 10f; // Duration of the movement

    private Vector3 startPosition;
    private Vector3 targetPosition;
    private float elapsedTime = 0f;

    void Start()
    {
        startPosition = transform.position;
        targetPosition = startPosition + new Vector3(-targetDistance, 0, 0);
    }

    void Update()
    {
        elapsedTime += Time.deltaTime;
        float t = Mathf.PingPong(elapsedTime / duration, 1);

        // Calculate the current position along the arc
        Vector3 currentPosition = Vector3.Lerp(startPosition, targetPosition, t);
        currentPosition.y += Mathf.Sin(t * Mathf.PI) * arcHeight;

        transform.position = currentPosition;
    }
}
