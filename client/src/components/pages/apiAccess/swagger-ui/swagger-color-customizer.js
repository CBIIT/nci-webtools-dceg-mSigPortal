// import React, { useEffect } from 'react';


// const SwaggerColorCustomizer = () => {
//   useEffect(() => {
//     const changeAllColors = () => {
//       // Select all elements with the style color rgb(211, 99, 99)
//       const elementsWithColor = document.querySelectorAll(
//         '[style*="color: rgb(211, 99, 99);"]'
//       );

//       let colorChanged = false;

//       elementsWithColor.forEach((span) => {
//         // Update the inline style to a new color
//         span.style.color = '#f2ee0f'; // color of the integer number
//         colorChanged = true;
//       });

//       // Stop the interval if colors have been changed
//       if (colorChanged) {
//         //console.log("Color change applied to integer values.");
//         clearInterval(intervalId);
//       }
//     };

//     // Set up an interval to retry color application every 1s
//     const intervalId = setInterval(() => {
//       changeAllColors();
//     }, 1000);

//     // Cleanup interval on component unmount
//     return () => clearInterval(intervalId);
//   }, []);

//   return null;
// };

// export default SwaggerColorCustomizer;

import React, { useEffect } from 'react';

const SwaggerColorCustomizer = () => {
  useEffect(() => {
    const highlightAllMatching = () => {
      // Get all code blocks inside tables (where response values usually appear)
      const elements = document.querySelectorAll('td code');

      elements.forEach((el) => {
        const style = window.getComputedStyle(el);
        const color = style.color;

        if (color === 'rgb(211, 99, 99)') {
          el.style.color = '#f2ee0f'; // yellow
          el.style.backgroundColor = 'black'; // optional: make it pop
          el.style.fontWeight = 'bold';
        }
      });
    };

    // Initial run
    highlightAllMatching();

    // Set up MutationObserver for dynamic changes
    const observer = new MutationObserver(() => {
      highlightAllMatching();
    });

    observer.observe(document.body, {
      childList: true,
      subtree: true,
    });

    // Clean up
    return () => observer.disconnect();
  }, []);

  return null;
};

export default SwaggerColorCustomizer;
