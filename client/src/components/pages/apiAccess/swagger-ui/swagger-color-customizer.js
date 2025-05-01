import React, { useEffect } from 'react';


const SwaggerColorCustomizer = () => {
  useEffect(() => {
    const changeAllColors = () => {
      // Select all elements with the style color rgb(211, 99, 99)
      const elementsWithColor = document.querySelectorAll(
        '[style*="color: rgb(211, 99, 99);"]'
      );

      let colorChanged = false;

      elementsWithColor.forEach((span) => {
        // Update the inline style to a new color
        span.style.color = '#f2ee0f'; // color of the integer number
        colorChanged = true;
      });
     
    };

    let tickCount = 0;
    const maxTicks = 300;

    const intervalId = setInterval(() => {
      changeAllColors();
      tickCount++;
      if (tickCount > maxTicks) clearInterval(intervalId);
    }, 1000);

    

    // Cleanup interval on component unmount
    return () => clearInterval(intervalId);
  }, []);

  return null;
};

export default SwaggerColorCustomizer;

