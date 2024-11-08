import React, { useEffect } from 'react';

const SwaggerScrollablePreEnhancer = () => {
    useEffect(() => {
      const addScrollingAccessibility = () => {
        // Select all <pre> elements that should be scrollable
        const preElements = document.querySelectorAll('pre[style*="overflow-x: auto"]');
  
        preElements.forEach((preElement) => {
          // Make the <pre> element focusable
          preElement.setAttribute('tabindex', '0');
  
          // Add keyboard scrolling functionality
          preElement.addEventListener('keydown', (event) => {
            const scrollAmount = 10; // Adjust for desired scroll sensitivity
  
            switch (event.key) {
              case 'ArrowUp':
                preElement.scrollTop -= scrollAmount;
                event.preventDefault();
                break;
              case 'ArrowDown':
                preElement.scrollTop += scrollAmount;
                event.preventDefault();
                break;
              case 'ArrowLeft':
                preElement.scrollLeft -= scrollAmount;
                event.preventDefault();
                break;
              case 'ArrowRight':
                preElement.scrollLeft += scrollAmount;
                event.preventDefault();
                break;
              default:
                break;
            }
          });
        });
      };
  
      // Retry to find and enhance <pre> elements every 500ms if not immediately available
      const intervalId = setInterval(() => {
        addScrollingAccessibility();
      }, 500);
  
      // Clear interval once all <pre> elements are enhanced
      if (document.querySelector('pre[style*="overflow-x: auto"]')) {
        clearInterval(intervalId);
      }
  
      // Cleanup interval on component unmount
      return () => clearInterval(intervalId);
    }, []);
  
    return null;
  };

export default SwaggerScrollablePreEnhancer;
