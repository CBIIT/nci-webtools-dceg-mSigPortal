// render `backtick` spans within a string as <code> elements
export function renderInlineCode(text) {
  if (typeof text !== 'string') return text;
  return text
    .split(/(`[^`]+`)/g)
    .map((part, i) =>
      part.startsWith('`') && part.endsWith('`') ? (
        <code key={i}>{part.slice(1, -1)}</code>
      ) : (
        part
      )
    );
}
